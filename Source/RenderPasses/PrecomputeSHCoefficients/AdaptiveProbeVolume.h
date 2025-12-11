#pragma once
#include "Falcor.h"
#include <vector>

#include "envMap_SH.h"

using namespace Falcor;

class AdaptiveProbeVolume : public Object
{
public:
    static ref<AdaptiveProbeVolume> create(ref<Device> pDevice);

    void startBuild(const ref<Scene>& pScene, float errorThreshold);
    bool hasPendingBatch() const { return !mPendingProbes.empty(); }
    void getPendingPositions(std::vector<float3> &position) const;

    // Updated: Input matches your structs directly
    void setProbeData(
        uint32_t batchIndex,
        const std::vector<float3>& coeffs,
        const std::vector<GradSHCoeff>& grads,
        const std::vector<HessianSHCoeff>& hessians
    );

    void finishBatch();
    void uploadToGPU();

    ref<Buffer> getNodeBuffer() const { return mpNodeBuffer; }
    ref<Buffer> getProbeBuffer() const { return mpProbeBuffer; }

    void printDebugInfo(const std::string& filename);

private:
    AdaptiveProbeVolume(ref<Device> pDevice);

    // ------------------------------------------------------------------
    // Internal Data Structures (No more flattening!)
    // ------------------------------------------------------------------
    struct ProbePoint
    {
        // Store RGB coefficients directly
        std::vector<float3> shCoeffs;

        // Store RGB gradients directly using your struct
        std::vector<GradSHCoeff> shGradients;

        float maxLambdaL2Norm = 0.0f;
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

    // ------------------------------------------------------------------
    // GPU Data Layouts
    // ------------------------------------------------------------------
    struct PackedNode
    {
        float3 minPoint;
        uint32_t childrenStartIndex;
        float3 maxPoint;
        uint32_t probeIndex;
    };

    // GPU-safe version of Gradient struct (using float4 for 16-byte alignment safety)
    struct GPUGradSHCoeff
    {
        float4 r; // xyz = gradient, w = pad
        float4 g;
        float4 b;
    };

    struct PackedProbe
    {
        // 9 Coefficients (RGB)
        // using float4 for alignment (xyz=rgb, w=pad)
        float4 coeffs[9];

        // 9 Gradients (RGB)
        GPUGradSHCoeff gradients[9];
    };

    // ... (Helpers and Members remain same) ...
    float computeHessianErrorNorm(const float* rawData); // Helper for error calc

    ref<Device> mpDevice;

    /*
     * DATA STRUCTURE RELATIONSHIP:
     * ----------------------------
     * 1. mPendingProbes (The Queue):
     * Stores a list of pairs <NodeIndex, ProbeIndex> representing the current batch
     * of work that just finished computing on the GPU.
     *
     * 2. mAdaptiveNodes (The Geometry):
     * Accessed via 'nodeIdx'. Represents the Octree box (min/max bounds).
     * We need this to calculate the voxel size (Delta x) for the error metric.
     *
     * 3. mProbePoints (The Data):
     * Accessed via 'node.probeIndex'. Represents the lighting samples.
     * We need this to retrieve the computed Hessian Error Norm (||lambda||).
     *
     * VISUALIZATION:
     * --------------
     * mPendingProbes (Queue)
     * +-----------+-----------+
     * idx: | Node ID: 5| Probe ID: 9|  <--- "pair"
     * +-----+-----+-----+-----+
     * |           |
     * v           v
     * mAdaptiveNodes      mProbePoints
     * +-------------+     +-----------------------+
     * 5 | Min: (0,0,0)|   9 | SH Coeffs: [0.5, ...]|
     * | Max: (1,1,1)|     | Hessians:  [...]     |
     * | ProbeIdx: 9 |---->| Error: 0.04          |
     * +-------------+     +-----------------------+
     */
    std::vector<AdaptiveNode> mAdaptiveNodes;
    std::vector<ProbePoint> mProbePoints;
    std::vector<std::pair<int, int>> mPendingProbes;

    float mCurrentThreshold = 0.01f;
    int mMaxLevel = 6;
    ref<Buffer> mpNodeBuffer;
    ref<Buffer> mpProbeBuffer;
};
