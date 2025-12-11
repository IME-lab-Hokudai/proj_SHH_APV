#include "AdaptiveProbeVolume.h"
#include <algorithm>
#include <queue>
#include <cmath>
#include <cstring> 
#include <fstream>
#include <iomanip>
// ------------------------------------------------------------------
// Helper: Analytic Eigenvalues for 3x3 Symmetric Matrix
// Used to compute curvature for the error metric (Eq. 22)
// ------------------------------------------------------------------
static void computeEigenvalues3x3(const float3x3& A, float& e1, float& e2, float& e3)
{
    const float ONE_THIRD = 1.0f / 3.0f;
    float m = (A[0][0] + A[1][1] + A[2][2]) * ONE_THIRD;

    // Shifted matrix K = A - mI
    float k00 = A[0][0] - m, k11 = A[1][1] - m, k22 = A[2][2] - m;
    float k01 = A[0][1], k02 = A[0][2], k12 = A[1][2];

    float q = 0.5f * (k00 * (k11 * k22 - k12 * k12) - k01 * (k01 * k22 - k12 * k02) + k02 * (k01 * k12 - k11 * k02));
    float p = (k00 * k00 + k11 * k11 + k22 * k22 + 2.0f * (k01 * k01 + k02 * k02 + k12 * k12)) / 6.0f;

    // Handle Identity-like matrices or numerical zero
    if (p < 1e-20f)
    {
        e1 = e2 = e3 = m;
        return;
    }

    float p_sqrt = std::sqrt(p);
    // Clamp for acos safety
    float det_val = q / (p * p_sqrt);
    det_val = std::max(-1.0f, std::min(1.0f, det_val));

    float phi = ONE_THIRD * std::acos(det_val);
    float two_sqrt_p = 2.0f * p_sqrt;
    float s = std::sin(phi), c = std::cos(phi);

    // Roots of the depressed cubic
    e1 = m + two_sqrt_p * c;
    e2 = m - two_sqrt_p * (c * 0.5f + s * 0.8660254f);
    e3 = m - two_sqrt_p * (c * 0.5f - s * 0.8660254f);
}

// ------------------------------------------------------------------
// Creation and Initialization
// ------------------------------------------------------------------

ref<AdaptiveProbeVolume> AdaptiveProbeVolume::create(ref<Device> pDevice)
{
    return ref<AdaptiveProbeVolume>(new AdaptiveProbeVolume(pDevice));
}

AdaptiveProbeVolume::AdaptiveProbeVolume(ref<Device> pDevice) : mpDevice(pDevice) {}

void AdaptiveProbeVolume::startBuild(const ref<Scene>& pScene, float errorThreshold, bool useRelativeError)
{
    // Reset all internal state
    mAdaptiveNodes.clear();
    mProbePoints.clear();
    mPendingProbes.clear();
    mCurrentThreshold = errorThreshold;
    mUseRelativeError = useRelativeError;

    // 1. Create Root Node (Covers entire scene AABB)
    // We treat the entire scene as a volume (no occupancy check needed)
    AdaptiveNode root;
    root.level = 0;
    root.minPoint = pScene->getSceneBounds().minPoint;
    root.maxPoint = pScene->getSceneBounds().maxPoint;

    // 2. Create the first placeholder probe at the root's center
    mProbePoints.emplace_back();
    root.probeIndex = 0;

    mAdaptiveNodes.push_back(root);

    // 3. Queue Root for processing (Level 0 Batch)
    mPendingProbes.push_back({0, 0});
}

// ------------------------------------------------------------------
// Build Pipeline: Step A (Prepare Data for GPU)
// ------------------------------------------------------------------

void AdaptiveProbeVolume::getPendingPositions(std::vector<float3> &positions) const
{
    positions.clear();
    if (mPendingProbes.empty())
        return;

    positions.reserve(mPendingProbes.size());

    // Extract center positions for all pending probes
    for (const auto& pair : mPendingProbes)
    {
        int nodeIdx = pair.first;
        float3 center = (mAdaptiveNodes[nodeIdx].minPoint + mAdaptiveNodes[nodeIdx].maxPoint) * 0.5f;
        positions.push_back(center);
    }
}

// ------------------------------------------------------------------
// Build Pipeline: Step B (Receive Data from CPU/GPU Math)
// ------------------------------------------------------------------

void AdaptiveProbeVolume::setProbeData(
    uint32_t batchIndex,
    const std::vector<float3>& coeffs,
    const std::vector<GradSHCoeff>& grads,
    const std::vector<HessianSHCoeff>& hessians
)
{
    // 1. Locate the correct probe in our internal list
    // 'batchIndex' matches the order in mPendingProbes/getPendingPositions
    int probeIdx = mPendingProbes[batchIndex].second;
    ProbePoint& pt = mProbePoints[probeIdx];

    // 2. Store SH Coefficients
    pt.shCoeffs = coeffs;
    float sumSqL = 0.0f;

    for (size_t i = 0; i < coeffs.size(); ++i)
    {
        // Accumulate squares for Vector Norm calculation
        sumSqL += math::dot(pt.shCoeffs[i], pt.shCoeffs[i]); // x^2 + y^2 + z^2
    }
    pt.coeffNorm = std::sqrt(sumSqL);
    // 3. Store Gradients (Copy directly)
    pt.shGradients = grads;

    // 4. Compute Hessian Error Metric (Eq. 22)
    // We treat R, G, and B as separate scalar fields and take the conservative maximum.
    float sumSquares = 0.0f;

    for (const auto& h : hessians)
    {
        // --- Red Channel ---
        float r1, r2, r3;
        computeEigenvalues3x3(h.r, r1, r2, r3);
        float maxR = std::max({std::abs(r1), std::abs(r2), std::abs(r3)});

        // --- Green Channel ---
        float g1, g2, g3;
        computeEigenvalues3x3(h.g, g1, g2, g3);
        float maxG = std::max({std::abs(g1), std::abs(g2), std::abs(g3)});

        // --- Blue Channel ---
        float b1, b2, b3;
        computeEigenvalues3x3(h.b, b1, b2, b3);
        float maxB = std::max({std::abs(b1), std::abs(b2), std::abs(b3)});

        // Conservative approach: Max curvature across any color channel dictates the error
        float maxLambda = std::max({maxR, maxG, maxB});

        // Accumulate for L2-Norm of the lambda vector
        sumSquares += maxLambda * maxLambda;
    }

    pt.maxLambdaL2Norm = std::sqrt(sumSquares);
}

// ------------------------------------------------------------------
// Build Pipeline: Step C (Subdivide and Schedule Next Batch)
// ------------------------------------------------------------------

void AdaptiveProbeVolume::finishBatch()
{
    std::vector<std::pair<int, int>> nextBatch;

    for (const auto& pair : mPendingProbes)
    {
        int nodeIdx = pair.first;

        // 1. SAFE READ: Access by index to check conditions
        // We do NOT store 'AdaptiveNode& node = ...' here because push_back below will invalidate it.
        if (mAdaptiveNodes[nodeIdx].level >= mMaxLevel)
            continue;

        ProbePoint& pt = mProbePoints[mAdaptiveNodes[nodeIdx].probeIndex];

        // E_abs calculation
        float3 diag = mAdaptiveNodes[nodeIdx].maxPoint - mAdaptiveNodes[nodeIdx].minPoint;
        float distSq = dot(diag, diag);
        float E_abs = (pt.maxLambdaL2Norm * distSq) * 0.5f;
        float finalError = E_abs;

        // 2. Apply Relative Error Logic if enabled
        if (mUseRelativeError)
        {
            // E_rel = E_abs / ||L||
            // Avoid division by zero for dark areas (use small epsilon)
            float L_norm = std::max(pt.coeffNorm, 1e-5f);
            finalError = E_abs / L_norm;
        }
        if (finalError > mCurrentThreshold)
        {
            // 2. CAPTURE DATA: Copy the values we need before modifying the vector.
            // This ensures 'minP' and 'size' stay valid even if the vector reallocates.
            float3 minP = mAdaptiveNodes[nodeIdx].minPoint;
            float3 maxP = mAdaptiveNodes[nodeIdx].maxPoint;
            float3 size = (maxP - minP) * 0.5f;
            int currentLevel = mAdaptiveNodes[nodeIdx].level;

            // 3. SAFE WRITE: Update parent flags using INDEX (stable)
            mAdaptiveNodes[nodeIdx].isLeaf = false;

            for (int k = 0; k < 8; ++k)
            {
                AdaptiveNode child;

                // Use the COPIED level variable
                child.level = currentLevel + 1;

                float3 offset = float3((k & 4) ? size.x : 0, (k & 2) ? size.y : 0, (k & 1) ? size.z : 0);
                child.minPoint = minP + offset;
                child.maxPoint = child.minPoint + size;

                // Create new probe
                mProbePoints.emplace_back();
                child.probeIndex = (uint32_t)mProbePoints.size() - 1;

                // 4. THE DANGER ZONE: This push_back might move the vector in memory
                mAdaptiveNodes.push_back(child);

                // 5. SAFE LINK: Re-access parent by INDEX to link the child
                // Even if the vector moved, 'nodeIdx' is still the correct offset.
                mAdaptiveNodes[nodeIdx].children[k] = (int)mAdaptiveNodes.size() - 1;

                nextBatch.push_back({(int)mAdaptiveNodes.size() - 1, (int)mProbePoints.size() - 1});
            }
        }
    }
    mPendingProbes = nextBatch;
}

// ------------------------------------------------------------------
// Finalization: Upload to GPU
// ------------------------------------------------------------------

void AdaptiveProbeVolume::uploadToGPU()
{
    //// 1. Pack Probes (Data Payload)
    //std::vector<PackedProbe> packedProbes(mProbePoints.size());

    //for (size_t i = 0; i < mProbePoints.size(); ++i)
    //{
    //    const ProbePoint& cpuPt = mProbePoints[i];
    //    PackedProbe& gpuPt = packedProbes[i];

    //    // Ensure we don't overflow fixed array (Order 3 = 9 coeffs)
    //    size_t count = std::min((size_t)9, cpuPt.shCoeffs.size());

    //    for (size_t k = 0; k < count; ++k)
    //    {
    //        // Pack Coeffs: float3 -> float4 (alignment safe)
    //        gpuPt.coeffs[k] = float4(cpuPt.shCoeffs[k], 0.0f);

    //        // Pack Gradients: float3 -> float4 (alignment safe)
    //        gpuPt.gradients[k].r = float4(cpuPt.shGradients[k].r, 0.0f);
    //        gpuPt.gradients[k].g = float4(cpuPt.shGradients[k].g, 0.0f);
    //        gpuPt.gradients[k].b = float4(cpuPt.shGradients[k].b, 0.0f);
    //    }
    //}

    //mpProbeBuffer = Buffer::createStructured(
    //    mpDevice,
    //    sizeof(PackedProbe),
    //    (uint32_t)packedProbes.size(),
    //    Resource::BindFlags::ShaderResource,
    //    Buffer::CpuAccess::None,
    //    packedProbes.data()
    //);

    //// 2. Pack Nodes (Breadth-First Flattening for cache locality)
    //std::vector<PackedNode> packedNodes;
    //packedNodes.reserve(mNodes.size()); // Pre-allocate to prevent reallocation during loop

    //std::vector<uint32_t> cpuToGpuIndex(mNodes.size());
    //std::queue<int> queue;
    //queue.push(0); // Start at Root

    //// Pass 1: Flatten nodes into GPU array
    //while (!queue.empty())
    //{
    //    int cpuIdx = queue.front();
    //    queue.pop();

    //    PackedNode gpuNode;
    //    gpuNode.minPoint = mNodes[cpuIdx].minPoint;
    //    gpuNode.maxPoint = mNodes[cpuIdx].maxPoint;
    //    gpuNode.probeIndex = mNodes[cpuIdx].probeIndex;
    //    gpuNode.childrenStartIndex = 0xFFFFFFFF; // Placeholder for leaf

    //    cpuToGpuIndex[cpuIdx] = (uint32_t)packedNodes.size();
    //    packedNodes.push_back(gpuNode);

    //    if (!mNodes[cpuIdx].isLeaf)
    //    {
    //        for (int k = 0; k < 8; ++k)
    //            queue.push(mNodes[cpuIdx].children[k]);
    //    }
    //}

    //// Pass 2: Patch children indices now that we know final positions
    //for (size_t i = 0; i < mNodes.size(); ++i)
    //{
    //    if (!mNodes[i].isLeaf)
    //    {
    //        uint32_t gpuParentIdx = cpuToGpuIndex[i];

    //        // Because we pushed children to the queue sequentially,
    //        // the GPU array will store them sequentially starting from child[0].
    //        int firstChildCpuIdx = mNodes[i].children[0];
    //        packedNodes[gpuParentIdx].childrenStartIndex = cpuToGpuIndex[firstChildCpuIdx];
    //    }
    //}

    //mpNodeBuffer = Buffer::createStructured(
    //    mpDevice,
    //    sizeof(PackedNode),
    //    (uint32_t)packedNodes.size(),
    //    Resource::BindFlags::ShaderResource,
    //    Buffer::CpuAccess::None,
    //    packedNodes.data()
    //);
}

void AdaptiveProbeVolume::printDebugInfo(const std::string& filename)
{
    std::ofstream out(filename);
    if (!out)
        return;

    out << "================================================================================\n";
    out << " ADAPTIVE PROBE VOLUME HIERARCHY\n";
    out << "================================================================================\n";
    out << " Nodes:     " << mAdaptiveNodes.size() << "\n";
    out << " Probes:    " << mProbePoints.size() << "\n";
    out << " Threshold: " << mCurrentThreshold << "\n";
    out << " Max Level: " << mMaxLevel << "\n";
    out << " Metric:    " << (mUseRelativeError ? "Relative (E_rel)" : "Absolute (E_abs)") << "\n";
    out << "================================================================================\n\n";

    std::function<void(int, std::string, bool)> printNode = [&](int nodeIdx, std::string prefix, bool isLast)
    {
        const AdaptiveNode& node = mAdaptiveNodes[nodeIdx];

        // 1. Calculate Geometry
        float3 diag = node.maxPoint - node.minPoint;
        float distSq = dot(diag, diag);
        float size = std::sqrt(distSq);

        // 2. Calculate Error (MATCHING finishBatch LOGIC)
        float error = 0.0f;
        float lambda = 0.0f;

        if (node.probeIndex < mProbePoints.size())
        {
            const ProbePoint& pt = mProbePoints[node.probeIndex];
            lambda = pt.maxLambdaL2Norm;

            // A. Absolute Error
            float E_abs = (lambda * distSq) * 0.5f;
            error = E_abs;

            // B. Relative Error (Must match finishBatch)
            if (mUseRelativeError)
            {
                float L_norm = std::max(pt.coeffNorm, 1e-5f);
                error = E_abs / L_norm;
            }
        }

        // 3. Print Structure
        out << prefix << (isLast ? "L-- " : "|-- ");
        out << "[Lvl " << node.level << "] ";

        // Format error string: "0.0050"
        std::stringstream ssErr;
        ssErr << std::fixed << std::setprecision(4) << error;
        std::string errStr = ssErr.str();

        if (!node.isLeaf)
        {
            // Case 1: SUBDIVIDED
            out << "SUBDIVIDED ";
            out << "(Err: " << errStr << " > " << mCurrentThreshold << ") ";
            out << "Size: " << std::setprecision(2) << size << "\n";
        }
        else
        {
            // Case 2: LEAF
            out << "LEAF       ";

            if (node.level == mMaxLevel)
            {
                out << "(Max Depth) ";
            }
            else if (error > mCurrentThreshold)
            {
                out << "BUG! (Err: " << errStr << " > Thr) ";
            }
            else
            {
                // explicit proof
                out << "(Err: " << errStr << " < " << mCurrentThreshold << ") ";
            }
            out << "Lambda: "  << std::setprecision(4) << lambda << "\n";
        }

        // 4. Recurse
        if (!node.isLeaf)
        {
            std::vector<int> childrenIdx;
            for (int i = 0; i < 8; ++i)
                if (node.children[i] != -1)
                    childrenIdx.push_back(node.children[i]);

            for (size_t i = 0; i < childrenIdx.size(); ++i)
            {
                std::string newPrefix = prefix + (isLast ? "    " : "|   ");
                printNode(childrenIdx[i], newPrefix, (i == childrenIdx.size() - 1));
            }
        }
    };

    if (!mAdaptiveNodes.empty())
    {
        printNode(0, "", true);
    }
    else
    {
        out << "Tree is empty.\n";
    }

    out.close();
}
