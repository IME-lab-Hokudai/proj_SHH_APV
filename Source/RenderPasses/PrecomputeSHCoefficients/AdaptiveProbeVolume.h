#pragma once
#include "Falcor.h"
#include <vector>

using namespace Falcor;

class AdaptiveProbeVolume : public Object
{
public:
    static ref<AdaptiveProbeVolume> create(ref<Device> pDevice);

    // ------------------------------------------------------------------
    // Phase 1: Initialization
    // ------------------------------------------------------------------
    // Resets the grid and creates the initial "Root Probe" to be computed.
    void startBuild(const ref<Scene>& pScene, float errorThreshold);

    // ------------------------------------------------------------------
    // Phase 2: The Build Loop (Called by RenderPass)
    // ------------------------------------------------------------------

    // Returns true if there is a batch of probes waiting for GPU computation.
    bool hasPendingBatch() const { return !mPendingProbes.empty(); }

    // Creates a buffer of float3 positions for the pending probes.
    // Bind this to your Compute Shader as "gProbePositions".
    void getPendingPositions(std::vector<float3> &positions);

    // Processes the results from the GPU.
    // 1. Reads back L/Grad/Hessian.
    // 2. Checks Error Metric.
    // 3. Subdivides nodes and generates the NEXT batch of probes for the next loop.
    // pResultBuffer struct: { float coeffs[9]; float3 grad[9]; float3x3 hess[9]; }
    void processBatchResults(const ref<Buffer>& pResultBuffer);

    // ------------------------------------------------------------------
    // Phase 3: Finalize
    // ------------------------------------------------------------------
    void uploadToGPU();

    ref<Buffer> getNodeBuffer() const { return mpNodeBuffer; }
    ref<Buffer> getProbeBuffer() const { return mpProbeBuffer; }

private:
    AdaptiveProbeVolume(ref<Device> pDevice);

    struct ProbePoint
    {
        std::vector<float> shCoeffs;
        std::vector<float3> shGradients;
        float hessianErrorNorm = 0.0f;
    };

    struct AdaptiveNode
    {
        bool isLeaf = true;
        int level = 0;
        float3 minPoint;
        float3 maxPoint;
        uint32_t probeIndex = 0;
        int children[8] = {-1};
    };

    // GPU Layouts
    struct PackedNode
    {
        float3 minPoint;
        uint32_t childrenStartIndex;
        float3 maxPoint;
        uint32_t probeIndex;
    };

    struct PackedProbe
    {
        float coeffs[9];
        float3 gradients[9];
    };

    // Helper
    float computeHessianErrorNorm(const float* rawData);

    ref<Device> mpDevice;

    std::vector<AdaptiveNode> mNodes;
    std::vector<ProbePoint> mProbePoints;

    // The Work Queue: Pairs of <NodeIndex, ProbeIndex> that need computing
    std::vector<std::pair<int, int>> mPendingProbes;

    float mCurrentThreshold = 0.01f;
    int mMaxLevel = 6;

    ref<Buffer> mpNodeBuffer;
    ref<Buffer> mpProbeBuffer;
};
