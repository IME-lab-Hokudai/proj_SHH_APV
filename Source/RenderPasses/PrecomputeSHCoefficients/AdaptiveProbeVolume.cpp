#include "AdaptiveProbeVolume.h"
#include <algorithm>
#include <queue>
#include <cmath>

// (Insert computeEigenvalues3x3 helper function here)
static void computeEigenvalues3x3(const float3x3& A, float& e1, float& e2, float& e3)
{
    const float ONE_THIRD = 1.0f / 3.0f;
    float m = (A[0][0] + A[1][1] + A[2][2]) * ONE_THIRD;

    float k00 = A[0][0] - m, k11 = A[1][1] - m, k22 = A[2][2] - m;
    float k01 = A[0][1], k02 = A[0][2], k12 = A[1][2];

    float q = 0.5f * (k00 * (k11 * k22 - k12 * k12) - k01 * (k01 * k22 - k12 * k02) + k02 * (k01 * k12 - k11 * k02));
    float p = (k00 * k00 + k11 * k11 + k22 * k22 + 2.0f * (k01 * k01 + k02 * k02 + k12 * k12)) / 6.0f;

    if (p < 1e-20f)
    {
        e1 = e2 = e3 = m;
        return;
    }

    float p_sqrt = std::sqrt(p);
    float det_val = std::max(-1.0f, std::min(1.0f, q / (p * p_sqrt)));

    float phi = ONE_THIRD * std::acos(det_val);
    float two_sqrt_p = 2.0f * p_sqrt;
    float s = std::sin(phi), c = std::cos(phi);

    e1 = m + two_sqrt_p * c;
    e2 = m - two_sqrt_p * (c * 0.5f + s * 0.8660254f);
    e3 = m - two_sqrt_p * (c * 0.5f - s * 0.8660254f);
}

ref<AdaptiveProbeVolume> AdaptiveProbeVolume::create(ref<Device> pDevice)
{
    return ref<AdaptiveProbeVolume>(new AdaptiveProbeVolume(pDevice));
}

AdaptiveProbeVolume::AdaptiveProbeVolume(ref<Device> pDevice) : mpDevice(pDevice) {}

void AdaptiveProbeVolume::startBuild(const ref<Scene>& pScene, float errorThreshold)
{
    mNodes.clear();
    mProbePoints.clear();
    mPendingProbes.clear();
    mCurrentThreshold = errorThreshold;

    // 1. Create Root Node
    AdaptiveNode root;
    root.level = 0;
    root.minPoint = pScene->getSceneBounds().minPoint;
    root.maxPoint = pScene->getSceneBounds().maxPoint;

    // 2. Create Placeholder Root Probe
    mProbePoints.emplace_back();
    root.probeIndex = 0;

    mNodes.push_back(root);

    // 3. Queue Root for processing (Level 0 Batch)
    mPendingProbes.push_back({0, 0});
}

void AdaptiveProbeVolume::getPendingPositions(std::vector<float3>& positions)
{
    positions.clear();
    positions.reserve(mPendingProbes.size());

    // Extract center positions for all pending probes
    for (const auto& pair : mPendingProbes)
    {
        int nodeIdx = pair.first;
        float3 center = (mNodes[nodeIdx].minPoint + mNodes[nodeIdx].maxPoint) * 0.5f;
        positions.push_back(center);
    }
}

void AdaptiveProbeVolume::processBatchResults(const ref<Buffer>& pResultBuffer)
{
    // Map GPU results to CPU
    // Stride = 117 floats (9 Coeffs + 27 Grads + 81 Hessians)
    //const float* pData = (const float*)pResultBuffer->map(Buffer::MapType::Read);
    const float* pData;

    std::vector<std::pair<int, int>> nextBatch;

    // Process every probe in the current batch
    for (size_t i = 0; i < mPendingProbes.size(); ++i)
    {
        int nodeIdx = mPendingProbes[i].first;
        int probeIdx = mPendingProbes[i].second;
        AdaptiveNode& node = mNodes[nodeIdx];
        ProbePoint& pt = mProbePoints[probeIdx];

        // --- 1. Store Data ---
        const float* pProbeData = pData + (i * 117);

        pt.shCoeffs.assign(pProbeData, pProbeData + 9);
        pt.shGradients.resize(9);
        memcpy(pt.shGradients.data(), pProbeData + 9, sizeof(float3) * 9);

        // Compute Error Metric from Hessian data (Offset 36)
        pt.hessianErrorNorm = computeHessianErrorNorm(pProbeData + 36);

        // --- 2. Check Subdivision ---
        if (node.level >= mMaxLevel)
            continue;

        // E_abs = ||lambda|| * ||Delta x||^2 / 2
        float3 diag = node.maxPoint - node.minPoint;
        float distSq = dot(diag, diag);
        float error = (pt.hessianErrorNorm * distSq) * 0.5f;

        if (error > mCurrentThreshold)
        {
            node.isLeaf = false;
            float3 size = (node.maxPoint - node.minPoint) * 0.5f;

            // Generate 8 children immediately
            for (int k = 0; k < 8; ++k)
            {
                AdaptiveNode child;
                child.level = node.level + 1;
                float3 offset = float3((k & 4) ? size.x : 0, (k & 2) ? size.y : 0, (k & 1) ? size.z : 0);
                child.minPoint = node.minPoint + offset;
                child.maxPoint = child.minPoint + size;

                // Create new uncomputed probe for child
                mProbePoints.emplace_back();
                child.probeIndex = (uint32_t)mProbePoints.size() - 1;

                mNodes.push_back(child);
                node.children[k] = (int)mNodes.size() - 1;

                // ** QUEUE FOR NEXT PASS **
                nextBatch.push_back({(int)mNodes.size() - 1, (int)mProbePoints.size() - 1});
            }
        }
    }

    pResultBuffer->unmap();

    // Replace current queue with the next level
    mPendingProbes = nextBatch;
}

float AdaptiveProbeVolume::computeHessianErrorNorm(const float* hessianData)
{
    float sumSquares = 0.0f;
    for (int k = 0; k < 9; ++k)
    {
        float3x3 H;
        memcpy(&H, hessianData + (k * 9), sizeof(float) * 9);

        float e1, e2, e3;
        computeEigenvalues3x3(H, e1, e2, e3);
        float maxLambda = std::max({std::abs(e1), std::abs(e2), std::abs(e3)});
        sumSquares += maxLambda * maxLambda;
    }
    return std::sqrt(sumSquares);
}

void AdaptiveProbeVolume::uploadToGPU()
{
    //// 1. Pack Probes
    //std::vector<PackedProbe> packedProbes(mProbePoints.size());
    //for (size_t i = 0; i < mProbePoints.size(); ++i)
    //{
    //    // Safe copy 9 coeffs
    //    size_t count = std::min((size_t)9, mProbePoints[i].shCoeffs.size());
    //    memcpy(packedProbes[i].coeffs, mProbePoints[i].shCoeffs.data(), sizeof(float) * count);
    //    memcpy(packedProbes[i].gradients, mProbePoints[i].shGradients.data(), sizeof(float3) * count);
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
    //std::vector<uint32_t> cpuToGpuIndex(mNodes.size());
    //std::queue<int> queue;
    //queue.push(0); // Root

    //while (!queue.empty())
    //{
    //    int cpuIdx = queue.front();
    //    queue.pop();

    //    PackedNode gpuNode;
    //    gpuNode.minPoint = mNodes[cpuIdx].minPoint;
    //    gpuNode.maxPoint = mNodes[cpuIdx].maxPoint;
    //    gpuNode.probeIndex = mNodes[cpuIdx].probeIndex;
    //    gpuNode.childrenStartIndex = 0xFFFFFFFF;

    //    cpuToGpuIndex[cpuIdx] = (uint32_t)packedNodes.size();
    //    packedNodes.push_back(gpuNode);

    //    if (!mNodes[cpuIdx].isLeaf)
    //    {
    //        for (int k = 0; k < 8; ++k)
    //            queue.push(mNodes[cpuIdx].children[k]);
    //    }
    //}

    //// Patch children indices
    //for (size_t i = 0; i < mNodes.size(); ++i)
    //{
    //    if (!mNodes[i].isLeaf)
    //    {
    //        uint32_t gpuParentIdx = cpuToGpuIndex[i];
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
