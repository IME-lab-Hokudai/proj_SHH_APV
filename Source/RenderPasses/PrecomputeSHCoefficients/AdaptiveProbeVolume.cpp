#include "AdaptiveProbeVolume.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>

// ------------------------------------------------------------------
// Helper: Analytic Eigenvalues for 3x3 Symmetric Matrix
// ------------------------------------------------------------------
static void computeEigenvalues3x3(const float3x3& A, float& e1, float& e2, float& e3)
{
    const float ONE_THIRD = 1.0f / 3.0f;
    float m = (A[0][0] + A[1][1] + A[2][2]) * ONE_THIRD;
    float k00 = A[0][0] - m, k11 = A[1][1] - m, k22 = A[2][2] - m;
    float k01 = A[0][1], k02 = A[0][2], k12 = A[1][2];
    float q = 0.5f * (k00 * (k11 * k22 - k12 * k12) - k01 * (k01 * k22 - k12 * k02) + k02 * (k01 * k12 - k11 * k02));
    float p = (k00 * k00 + k11 * k11 + k22 * k22 + 2.0f * (k01 * k01 + k02 * k02 + k12 * k12)) / 6.0f;

    if (p < 1e-20f) { e1 = e2 = e3 = m; return; }

    float p_sqrt = std::sqrt(p);
    float det_val = q / (p * p_sqrt);
    det_val = std::max(-1.0f, std::min(1.0f, det_val));
    float phi = ONE_THIRD * std::acos(det_val);
    float two_sqrt_p = 2.0f * p_sqrt;
    float s = std::sin(phi), c = std::cos(phi);

    e1 = m + two_sqrt_p * c;
    e2 = m - two_sqrt_p * (c * 0.5f + s * 0.8660254f);
    e3 = m - two_sqrt_p * (c * 0.5f - s * 0.8660254f);
}

// ------------------------------------------------------------------
// Lifecycle
// ------------------------------------------------------------------

ref<AdaptiveProbeVolume> AdaptiveProbeVolume::create(ref<Device> pDevice)
{
    return ref<AdaptiveProbeVolume>(new AdaptiveProbeVolume(pDevice));
}

AdaptiveProbeVolume::AdaptiveProbeVolume(ref<Device> pDevice) : mpDevice(pDevice) {}

// ------------------------------------------------------------------
// Build Logic
// ------------------------------------------------------------------

void AdaptiveProbeVolume::startBuild(const ref<Scene>& pScene, float errorThreshold, bool useRelativeError)
{
    mProbes.clear();
    mCorners.clear();
    mPendingNewCorners.clear();
    mProbesPendingCheck.clear();

    mCurrentThreshold = errorThreshold;
    mUseRelativeError = useRelativeError;

    auto bounds = pScene->getSceneBounds();

    // 1. Create Root Probe
    Probe root;
    root.level = 0;
    root.minPoint = bounds.minPoint;
    root.maxPoint = bounds.maxPoint;

    float3 size = root.maxPoint - root.minPoint;

    // 2. Create the initial 8 Corners
    for (int i = 0; i < 8; ++i)
    {
        float3 offset = float3((i & 4) ? size.x : 0, (i & 2) ? size.y : 0, (i & 1) ? size.z : 0);

        Corner c;
        c.position = root.minPoint + offset;

        mCorners.push_back(c);
        int idx = (int)mCorners.size() - 1;

        root.corners[i] = idx;
        mPendingNewCorners.push_back(idx);
    }

    mProbes.push_back(root);
    mProbesPendingCheck.push_back(0);
}

void AdaptiveProbeVolume::getPendingPositions(std::vector<float3>& positions) const
{
    positions.clear();
    positions.reserve(mPendingNewCorners.size());
    for (int idx : mPendingNewCorners)
    {
        positions.push_back(mCorners[idx].position);
    }
}

void AdaptiveProbeVolume::setCornerData(
    uint32_t batchIndex,
    const std::vector<float3>& coeffs,
    const std::vector<GradSHCoeff>& grads,
    const std::vector<float3x3>& hessians // Now receives Luminance Hessians
)
{
    // 1. Locate Probe
    if (batchIndex >= mPendingNewCorners.size()) return;

    int cornerIdx = mPendingNewCorners[batchIndex];
    Corner& c = mCorners[cornerIdx];

    // -----------------------------------------------------------
    // PART A: RUNTIME DATA (Store RGB)
    // -----------------------------------------------------------
    c.shCoeffs = coeffs;       // Full RGB for color rendering
    c.shGradients = grads;     // Full RGB for color interpolation

    // -----------------------------------------------------------
    // PART B: CONSTRUCTION METRICS (Use Luminance)
    // -----------------------------------------------------------
    if (mUseRelativeError) {
        const float3 kLuma = float3(0.2126f, 0.7152f, 0.0722f);

        // 1. Calculate Norm of LUMINANCE Coefficients (Vector L)
        // ||L|| = sqrt( sum( (c_i . luma)^2 ) )
        float sumSqL = 0.0f;
        for (const auto& c : coeffs)
        {
            // Project this band's RGB coefficient to Luminance
            float lum = dot(c, kLuma);
            sumSqL += lum * lum;
        }
        c.coeffVecL2 = std::sqrt(sumSqL);
    }

    // 2. Calculate Curvature of LUMINANCE Field (Hessian)
    float sumSquares = 0.0f;
    for (const auto& h : hessians)
    {
        // Compute Eigenvalues of the Luminance Hessian
        float e1, e2, e3;
        computeEigenvalues3x3(h, e1, e2, e3);

        // Max curvature
        float maxLambda = std::max({ std::abs(e1), std::abs(e2), std::abs(e3) });
        sumSquares += maxLambda * maxLambda;
    }

    c.maxLambdaVecL2 = std::sqrt(sumSquares);
}

void AdaptiveProbeVolume::finishBatch()
{
    // Clear pending list; we will fill it with the NEXT generation's corners
    mPendingNewCorners.clear();

    std::vector<int> nextProbesToCheck;

    for (int probeIdx : mProbesPendingCheck)
    {
        if (mProbes[probeIdx].level >= mMaxLevel) continue;

        // --------------------------------------------------
        // Step 4: Calculate Error for this Probe
        // --------------------------------------------------
        float3 diag = mProbes[probeIdx].maxPoint - mProbes[probeIdx].minPoint;
        float distSq = dot(diag, diag);

        float totalError = 0.0f;

        for (int k = 0; k < 8; ++k)
        {
            int cIdx = mProbes[probeIdx].corners[k];
            const Corner& c = mCorners[cIdx];

            // E_abs = 0.5 * Lambda * dx^2
            float e = (c.maxLambdaVecL2 * distSq) * 0.5f;

            if (mUseRelativeError)
            {
                float norm = std::max(c.coeffVecL2, 1e-5f);
                totalError += (e / norm);
            }
            else
            {
                totalError += e;
            }
        }
        float avgProbeError = totalError / 8.0f;

        // --------------------------------------------------
        // Step 5: Subdivide
        // --------------------------------------------------
        if (avgProbeError > mCurrentThreshold)
        {
            mProbes[probeIdx].isLeaf = false;

            float3 minP = mProbes[probeIdx].minPoint;
            float3 maxP = mProbes[probeIdx].maxPoint;
            float3 center = (minP + maxP) * 0.5f;

            // Helper: Create new corner and add to pending list
            auto addNewCorner = [&](float3 pos) -> int {
                Corner c; c.position = pos;
                mCorners.push_back(c);
                int idx = (int)mCorners.size() - 1;
                mPendingNewCorners.push_back(idx);
                return idx;
                };

            // A. Generate 19 NEW Corners
            int c_center = addNewCorner(center);

            int c_fX0 = addNewCorner({ minP.x, center.y, center.z }); int c_fX1 = addNewCorner({ maxP.x, center.y, center.z });
            int c_fY0 = addNewCorner({ center.x, minP.y, center.z }); int c_fY1 = addNewCorner({ center.x, maxP.y, center.z });
            int c_fZ0 = addNewCorner({ center.x, center.y, minP.z }); int c_fZ1 = addNewCorner({ center.x, center.y, maxP.z });

            int c_eX_Y0Z0 = addNewCorner({ center.x, minP.y, minP.z }); int c_eX_Y1Z0 = addNewCorner({ center.x, maxP.y, minP.z });
            int c_eX_Y0Z1 = addNewCorner({ center.x, minP.y, maxP.z }); int c_eX_Y1Z1 = addNewCorner({ center.x, maxP.y, maxP.z });

            int c_eY_X0Z0 = addNewCorner({ minP.x, center.y, minP.z }); int c_eY_X1Z0 = addNewCorner({ maxP.x, center.y, minP.z });
            int c_eY_X0Z1 = addNewCorner({ minP.x, center.y, maxP.z }); int c_eY_X1Z1 = addNewCorner({ maxP.x, center.y, maxP.z });

            int c_eZ_X0Y0 = addNewCorner({ minP.x, minP.y, center.z }); int c_eZ_X1Y0 = addNewCorner({ maxP.x, minP.y, center.z });
            int c_eZ_X0Y1 = addNewCorner({ minP.x, maxP.y, center.z }); int c_eZ_X1Y1 = addNewCorner({ maxP.x, maxP.y, center.z });

            // B. Retrieve 8 OLD Corners
            const int* P = mProbes[probeIdx].corners;

            // C. Map to Grid
            int grid[3][3][3];
            grid[0][0][0] = P[0]; grid[0][0][2] = P[1]; grid[0][2][0] = P[2]; grid[0][2][2] = P[3];
            grid[2][0][0] = P[4]; grid[2][0][2] = P[5]; grid[2][2][0] = P[6]; grid[2][2][2] = P[7];
            grid[1][1][1] = c_center;
            grid[0][1][1] = c_fX0; grid[2][1][1] = c_fX1;
            grid[1][0][1] = c_fY0; grid[1][2][1] = c_fY1;
            grid[1][1][0] = c_fZ0; grid[1][1][2] = c_fZ1;
            grid[1][0][0] = c_eX_Y0Z0; grid[1][2][0] = c_eX_Y1Z0; grid[1][0][2] = c_eX_Y0Z1; grid[1][2][2] = c_eX_Y1Z1;
            grid[0][1][0] = c_eY_X0Z0; grid[2][1][0] = c_eY_X1Z0; grid[0][1][2] = c_eY_X0Z1; grid[2][1][2] = c_eY_X1Z1;
            grid[0][0][1] = c_eZ_X0Y0; grid[2][0][1] = c_eZ_X1Y0; grid[0][2][1] = c_eZ_X0Y1; grid[2][2][1] = c_eZ_X1Y1;

            // D. Create 8 Child Probes
            for (int k = 0; k < 8; ++k)
            {
                Probe child;
                child.level = mProbes[probeIdx].level + 1;

                float3 childSize = (maxP - minP) * 0.5f;
                float3 offset = float3((k & 4) ? childSize.x : 0, (k & 2) ? childSize.y : 0, (k & 1) ? childSize.z : 0);
                child.minPoint = minP + offset;
                child.maxPoint = child.minPoint + childSize;

                int startX = (k & 4) ? 1 : 0;
                int startY = (k & 2) ? 1 : 0;
                int startZ = (k & 1) ? 1 : 0;

                for (int c = 0; c < 8; ++c)
                {
                    int dx = (c & 4) ? 1 : 0;
                    int dy = (c & 2) ? 1 : 0;
                    int dz = (c & 1) ? 1 : 0;
                    child.corners[c] = grid[startX + dx][startY + dy][startZ + dz];
                }

                mProbes.push_back(child);
                mProbes[probeIdx].children[k] = (int)mProbes.size() - 1;
                nextProbesToCheck.push_back((int)mProbes.size() - 1);
            }
        }
    }
    mProbesPendingCheck = nextProbesToCheck;
}

// ------------------------------------------------------------------
// CPU Traversal Logic (Matches Shader Logic)
// ------------------------------------------------------------------
int AdaptiveProbeVolume::traverseOctreeCPU(float3 pos) const
{
    if (mProbes.empty()) return -1;

    // 1. Bounds Check: If position is outside the volume, return -1
    const Probe& root = mProbes[0];
    if (pos.x < root.minPoint.x || pos.y < root.minPoint.y || pos.z < root.minPoint.z ||
        pos.x > root.maxPoint.x || pos.y > root.maxPoint.y || pos.z > root.maxPoint.z)
    {
        return -1;
    }

    int currentIdx = 0; // Start at Root

    // Safety loop to prevent infinite recursion
    for (int i = 0; i < mMaxLevel + 5; ++i)
    {
        const Probe& p = mProbes[currentIdx];

        // Found a Leaf
        if (p.children[0] == -1)
        {
            return currentIdx;
        }

        // Determine Octant
        float3 center = (p.minPoint + p.maxPoint) * 0.5f;
        int octant = 0;
        if (pos.x >= center.x) octant |= 4; // Bit 2 = X
        if (pos.y >= center.y) octant |= 2; // Bit 1 = Y
        if (pos.z >= center.z) octant |= 1; // Bit 0 = Z

        int childIdx = p.children[octant];

        // If tree is incomplete, stop here
        if (childIdx == -1) return currentIdx;

        currentIdx = childIdx;
    }
    return currentIdx;
}

// ------------------------------------------------------------------
// Neighbor Computation (Fixes T-Junctions)
// ------------------------------------------------------------------
void AdaptiveProbeVolume::computeNeighbors()
{
    // Directions corresponding to the 6 faces: -X, +X, -Y, +Y, -Z, +Z
    const float3 DIRECTIONS[6] = {
        float3(-1, 0, 0), float3(1, 0, 0),
        float3(0, -1, 0), float3(0, 1, 0),
        float3(0, 0, -1), float3(0, 0, 1)
    };

    for (size_t i = 0; i < mProbes.size(); ++i)
    {
        Probe& p = mProbes[i];

        // Initialize neighbors to -1 (No neighbor / No blend needed)
        for (int k = 0; k < 6; ++k) p.coarseNeighbors[k] = -1;

        // We only care about LEAF nodes
        if (p.children[0] != -1) continue;

        float3 size = p.maxPoint - p.minPoint;
        float3 center = (p.minPoint + p.maxPoint) * 0.5f;
        int totalCoarseNeighbors = 0;
        // Check all 6 faces
        for (int face = 0; face < 6; ++face)
        {
            // Step 5% of the size past the wall to ensure we enter the neighbor's voxel
            float3 dir = DIRECTIONS[face];
            float3 wallPos = center + dir * (size * 0.5f);
            float3 nudge = dir * (size * 0.1f); // Push 10% outwards
            float3 checkPos = wallPos + nudge;

            int neighborIdx = traverseOctreeCPU(checkPos);

            // If we found a valid neighbor (that is not ourselves)
            if (neighborIdx != -1 && neighborIdx != (int)i)
            {
                const Probe& neighbor = mProbes[neighborIdx];
                float3 neighborSize = neighbor.maxPoint - neighbor.minPoint;

                // THE RULE: Only store if neighbor is LARGER (Coarser)
                // Use 1.1x tolerance to avoid floating point equality issues
                if (neighborSize.x > size.x * 1.1f)
                {
                    p.coarseNeighbors[face] = neighborIdx;
                    totalCoarseNeighbors++;
                }
            }
        }
        std::cout << "Found " << totalCoarseNeighbors << " coarse neighbor links." << std::endl;
    }
}

void AdaptiveProbeVolume::uploadToGPU()
{
    // Ensure neighbors are computed before uploading!
    computeNeighbors();

    // 1. Pack Probes (Tree Topology)
    std::vector<GPUProbe> gpuProbes;
    gpuProbes.reserve(mProbes.size());

    for (const auto& p : mProbes)
    {
        GPUProbe gp;
        gp.minPoint = p.minPoint;
        gp.maxPoint = p.maxPoint;
        gp.pad1 = 0; gp.pad2 = 0;

        for (int i = 0; i < 8; ++i) gp.children[i] = p.children[i];
        for (int i = 0; i < 8; ++i) gp.corners[i] = p.corners[i];

        // --------------------------------------------------------
        // NEW: Copy Coarse Neighbor Indices
        // --------------------------------------------------------
        for (int k = 0; k < 6; ++k) gp.coarseNeighbors[k] = p.coarseNeighbors[k];
        gp.pad3[0] = 0; // Padding maintenance
        gp.pad3[1] = 0;

        gpuProbes.push_back(gp);
    }

    mpProbeBuffer = mpDevice->createStructuredBuffer(
        sizeof(GPUProbe),
        (uint32_t)gpuProbes.size(),
        ResourceBindFlags::ShaderResource,
        MemoryType::DeviceLocal,
        gpuProbes.data()
    );

    // 2. Pack Corners (Physics Data)
    std::vector<GPUCorner> gpuCorners;
    gpuCorners.reserve(mCorners.size());

    for (const auto& c : mCorners)
    {
        GPUCorner gc;

        // Safety: ensure we don't read out of bounds if corner has fewer bands
        int numBands = std::min((int)c.shCoeffs.size(), 9);

        for (int i = 0; i < 9; ++i)
        {
            if (i < numBands)
            {
                // Value
                gc.coeffs[i] = float4(c.shCoeffs[i], 0.0f);

                // Gradients
                gc.gradR[i] = float4(c.shGradients[i].r, 0.0f);
                gc.gradG[i] = float4(c.shGradients[i].g, 0.0f);
                gc.gradB[i] = float4(c.shGradients[i].b, 0.0f);
            }
            else
            {
                gc.coeffs[i] = float4(0);
                gc.gradR[i] = float4(0);
                gc.gradG[i] = float4(0);
                gc.gradB[i] = float4(0);
            }
        }
        gpuCorners.push_back(gc);
    }

    mpCornerBuffer = mpDevice->createStructuredBuffer(
        sizeof(GPUCorner),
        (uint32_t)gpuCorners.size(),
        ResourceBindFlags::ShaderResource,
        MemoryType::DeviceLocal,
        gpuCorners.data()
    );
}

void AdaptiveProbeVolume::printDebugInfo(const std::string& filename)
{
    std::ofstream out(filename);
    if (!out) return;

    out << "================================================================================\n";
    out << " ADAPTIVE PROBE VOLUME HIERARCHY (Vertex-Sharing)\n";
    out << "================================================================================\n";
    out << " Probes:    " << mProbes.size() << "\n";
    out << " Corners:   " << mCorners.size() << "\n";
    out << " Threshold: " << mCurrentThreshold << "\n";
    out << " Max Level: " << mMaxLevel << "\n";
    out << " Metric:    " << (mUseRelativeError ? "Relative (E_rel)" : "Absolute (E_abs)") << "\n";
    out << "================================================================================\n\n";

    // Recursive Lambda for Tree Traversal
    std::function<void(int, std::string, bool)> printProbe = [&](int probeIdx, std::string prefix, bool isLast)
        {
            const Probe& probe = mProbes[probeIdx];

            // 1. Re-Calculate Geometry
            float3 diag = probe.maxPoint - probe.minPoint;
            float distSq = dot(diag, diag);
            float size = std::sqrt(distSq);

            // 2. Re-Calculate Error (Must match finishBatch logic exactly)
            float totalError = 0.0f;
            float maxLambdaInProbe = 0.0f;

            for (int k = 0; k < 8; ++k)
            {
                int cIdx = probe.corners[k];
                const Corner& c = mCorners[cIdx];

                // Track max curvature in this voxel for display
                maxLambdaInProbe = std::max(maxLambdaInProbe, c.maxLambdaVecL2);

                // Re-compute error
                float e = (c.maxLambdaVecL2 * distSq) * 0.5f;

                if (mUseRelativeError)
                {
                    float norm = std::max(c.coeffVecL2, 1e-5f);
                    totalError += (e / norm);
                }
                else
                {
                    totalError += e;
                }
            }
            float avgError = totalError / 8.0f;

            // 3. Format Output
            out << prefix << (isLast ? "L-- " : "|-- ");
            out << "[Lvl " << probe.level << "] ";

            // Format error: "0.0050"
            std::stringstream ssErr;
            ssErr << std::fixed << std::setprecision(4) << avgError;
            std::string errStr = ssErr.str();

            if (!probe.isLeaf)
            {
                // Case 1: SUBDIVIDED
                out << "SUBDIVIDED ";
                out << "(Err: " << errStr << " > " << mCurrentThreshold << ") ";
                out << "Size: " << std::setprecision(2) << size << " ";
                out << "MaxLambda: " << std::setprecision(3) << maxLambdaInProbe << "\n";

                // Recurse children
                std::vector<int> childrenIdx;
                for (int i = 0; i < 8; ++i)
                    if (probe.children[i] != -1)
                        childrenIdx.push_back(probe.children[i]);

                for (size_t i = 0; i < childrenIdx.size(); ++i)
                {
                    std::string newPrefix = prefix + (isLast ? "    " : "|   ");
                    printProbe(childrenIdx[i], newPrefix, (i == childrenIdx.size() - 1));
                }
            }
            else
            {
                // Case 2: LEAF
                out << "LEAF       ";

                if (probe.level == mMaxLevel)
                {
                    out << "(Max Depth) ";
                }
                else if (avgError > mCurrentThreshold)
                {
                    // This happens if finishBatch hasn't processed this leaf yet, or logic bug
                    out << "PENDING? (Err: " << errStr << " > Thr) ";
                }
                else
                {
                    out << "(Err: " << errStr << " < " << mCurrentThreshold << ") ";
                }
                out << "MaxLambda: " << std::setprecision(3) << maxLambdaInProbe << "\n";
            }
        };

    if (!mProbes.empty())
    {
        printProbe(0, "", true);
    }
    else
    {
        out << "Tree is empty.\n";
    }

    out.close();
}

void AdaptiveProbeVolume::saveToFile(const std::string& filename) const
{
    std::ofstream out(filename);
    if (!out)
    {
        logError("Failed to open file for saving: " + filename);
        return;
    }

    out << std::fixed << std::setprecision(8);

    // 1. Header & Config
    out << "ADAPTIVE_GRID_V2\n";
    out << mCurrentThreshold << " " << mMaxLevel << " " << (mUseRelativeError ? 1 : 0) << "\n";

    // 2. Corners
    out << "NUM_CORNERS " << mCorners.size() << "\n";
    for (const auto& c : mCorners)
    {
        // Basic Data
        out << c.position.x << " " << c.position.y << " " << c.position.z << " ";
        out << c.maxLambdaVecL2 << " " << c.coeffVecL2 << "\n";

        // SH Coeffs
        out << c.shCoeffs.size() << "\n";
        for (const auto& val : c.shCoeffs)
        {
            out << val.x << " " << val.y << " " << val.z << " ";
        }
        out << "\n";

        // SH Gradients
        out << c.shGradients.size() << "\n";
        for (const auto& g : c.shGradients)
        {
            // Assuming GradSHCoeff has r, g, b (float3)
            out << g.r.x << " " << g.r.y << " " << g.r.z << " ";
            out << g.g.x << " " << g.g.y << " " << g.g.z << " ";
            out << g.b.x << " " << g.b.y << " " << g.b.z << " ";
        }
        out << "\n";
    }

    // 3. Probes
    out << "NUM_PROBES " << mProbes.size() << "\n";
    for (const auto& p : mProbes)
    {
        out << (p.isLeaf ? 1 : 0) << " " << p.level << "\n";

        // Bounds
        out << p.minPoint.x << " " << p.minPoint.y << " " << p.minPoint.z << " ";
        out << p.maxPoint.x << " " << p.maxPoint.y << " " << p.maxPoint.z << "\n";

        // Corners
        for (int i = 0; i < 8; ++i) out << p.corners[i] << " ";
        out << "\n";

        // Children
        for (int i = 0; i < 8; ++i) out << p.children[i] << " ";
        out << "\n";
    }

    out.close();
    logInfo("Successfully saved AdaptiveProbeVolume to " + filename);
}

void AdaptiveProbeVolume::loadFromFile(const std::string& filename)
{
    std::ifstream in(filename);
    if (!in)
    {
        logError("Failed to open file for loading: " + filename);
        return;
    }

    // Clear existing
    mProbes.clear();
    mCorners.clear();
    mPendingNewCorners.clear();
    mProbesPendingCheck.clear();

    std::string header;
    in >> header;
    if (header != "ADAPTIVE_GRID_V2")
    {
        logError("Invalid file format or version mismatch: " + filename);
        return;
    }

    int useRelErrInt;
    in >> mCurrentThreshold >> mMaxLevel >> useRelErrInt;
    mUseRelativeError = (useRelErrInt != 0);

    // 2. Load Corners
    size_t numCorners;
    std::string tag;
    in >> tag >> numCorners;

    mCorners.resize(numCorners);
    for (size_t i = 0; i < numCorners; ++i)
    {
        Corner& c = mCorners[i];
        in >> c.position.x >> c.position.y >> c.position.z;
        in >> c.maxLambdaVecL2 >> c.coeffVecL2;

        // SH Coeffs
        size_t numCoeffs;
        in >> numCoeffs;
        c.shCoeffs.resize(numCoeffs);
        for (size_t k = 0; k < numCoeffs; ++k)
        {
            in >> c.shCoeffs[k].x >> c.shCoeffs[k].y >> c.shCoeffs[k].z;
        }

        // SH Gradients
        size_t numGrads;
        in >> numGrads;
        c.shGradients.resize(numGrads);
        for (size_t k = 0; k < numGrads; ++k)
        {
            in >> c.shGradients[k].r.x >> c.shGradients[k].r.y >> c.shGradients[k].r.z;
            in >> c.shGradients[k].g.x >> c.shGradients[k].g.y >> c.shGradients[k].g.z;
            in >> c.shGradients[k].b.x >> c.shGradients[k].b.y >> c.shGradients[k].b.z;
        }
    }

    // 3. Load Probes
    size_t numProbes;
    in >> tag >> numProbes;

    mProbes.resize(numProbes);
    for (size_t i = 0; i < numProbes; ++i)
    {
        Probe& p = mProbes[i];
        int isLeafInt;
        in >> isLeafInt >> p.level;
        p.isLeaf = (isLeafInt != 0);

        in >> p.minPoint.x >> p.minPoint.y >> p.minPoint.z;
        in >> p.maxPoint.x >> p.maxPoint.y >> p.maxPoint.z;

        for (int k = 0; k < 8; ++k) in >> p.corners[k];
        for (int k = 0; k < 8; ++k) in >> p.children[k];
    }

    in.close();
    logInfo("Successfully loaded AdaptiveProbeVolume from " + filename);
}
