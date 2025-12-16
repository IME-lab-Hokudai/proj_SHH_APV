#pragma once
#include "Falcor.h"
#include <vector>
#include "envMap_SH.h"

using namespace Falcor;

class AdaptiveProbeVolume : public Object
{
public:
    static ref<AdaptiveProbeVolume> create(ref<Device> pDevice);

    void startBuild(const ref<Scene>& pScene, float errorThreshold, bool useRelativeError = false);

    // ----------------------------------------------------------------
    // Batch Processing Interface
    // ----------------------------------------------------------------

    // Check if there are new corners waiting for SH/Grad calculation
    bool hasPendingBatch() const { return !mPendingNewCorners.empty(); }

    // Get world positions of strictly the NEW corners that need calculation
    void getPendingPositions(std::vector<float3>& positions) const;

    // Fill physics data for a SPECIFIC corner in the current batch.
    // batchIndex: The index in mPendingNewCorners to update.
    // coeffs: All SH coefficients for this corner.
    // grads: All gradients for this corner.
    void setCornerData(
        uint32_t batchIndex,
        const std::vector<float3>& coeffs,
        const std::vector<GradSHCoeff>& grads,
        const std::vector<HessianSHCoeff>& hessians
    );

    // Calculate Error per Probe, Subdivide if needed, Generate new Corners
    void finishBatch();

    // ----------------------------------------------------------------
    // Resources & Debug
    // ----------------------------------------------------------------
    void uploadToGPU();
    void printDebugInfo(const std::string& filename);

    ref<Buffer> getNodeBuffer() const { return mpProbeBuffer; }
    ref<Buffer> getProbeBuffer() const { return mpCornerBuffer; }

private:
    AdaptiveProbeVolume(ref<Device> pDevice);

    // ----------------------------------------------------------------
    // Internal Data Structures
    // ----------------------------------------------------------------

    struct Corner
    {
        float3 position;

        // Physics Data (Full bands)
        std::vector<float3> shCoeffs;
        std::vector<GradSHCoeff> shGradients;

        float maxLambda = 0.0f; // Curvature
        float coeffNorm = 0.0f; // ||L||
    };

    struct Probe
    {
        bool isLeaf = true;
        int level = 0;

        float3 minPoint;
        float3 maxPoint;

        // Indices into mCorners (0..7)
        int corners[8] = { -1 };

        // Indices into mProbes (Children)
        int children[8] = { -1 };
    };

    ref<Device> mpDevice;

    std::vector<Probe> mProbes;
    std::vector<Corner> mCorners;

    // QUEUES
    std::vector<int> mPendingNewCorners;
    std::vector<int> mProbesPendingCheck;

    // Settings
    float mCurrentThreshold = 0.01f;
    int mMaxLevel = 5;
    bool mUseRelativeError = false;

    ref<Buffer> mpProbeBuffer;
    ref<Buffer> mpCornerBuffer;
};
