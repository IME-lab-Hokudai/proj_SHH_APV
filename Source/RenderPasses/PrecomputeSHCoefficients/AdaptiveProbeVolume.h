#pragma once
#include "Falcor.h"
#include <vector>
#include "envMap_SH.h"

using namespace Falcor;

class AdaptiveProbeVolume : public Object
{
public:
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

        int coarseNeighbors[6];
    };
    struct Corner
    {
        float3 position;

        // Physics Data (Full bands)
        std::vector<float3> shCoeffs;
        std::vector<GradSHCoeff> shGradients;

        float maxLambdaVecL2 = 0.0f; // Curvature
        float coeffVecL2 = 0.0f; // ||L||
        bool isValid;
        std::vector<float> distMean;   // New
        std::vector<float> distMeanSq; // New
        float constraintWeight = 0.0f;
    };

    //
// Add this helper to quantize positions (merges points closer than 0.1mm)
    struct CornerKey {
        int x, y, z;

        bool operator==(const CornerKey& other) const {
            return x == other.x && y == other.y && z == other.z;
        }
    };

    struct CornerKeyHasher {
        std::size_t operator()(const CornerKey& k) const {
            // Simple hash combination
            return ((std::hash<int>()(k.x) ^ (std::hash<int>()(k.y) << 1)) >> 1) ^ (std::hash<int>()(k.z) << 1);
        }
    };

    void interpolateRobust_CPU(int probeIdx, float3 pos, std::vector<float3>& outCoeffs);

    void fetchSmartInterpolatedSH(int probeIdx, float3 pos, float3* finalCoeffs);

    void interpolateHermite_CPU(int coarseProbeIdx, float3 pos, std::vector<float3>& outCoeffs, std::vector<GradSHCoeff>& outGrads);

    void constrainHangingNodes();
    void constrainHangingNodesHermite();

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
        const std::vector<float3x3>& hessians,
        const std::vector<float>& distMeans,
        const std::vector<float>& distMeanSqs
    );

    // Calculate Error per Probe, Subdivide if needed, Generate new Corners
    void finishBatch();

    int traverseOctreeCPU(float3 pos) const;

    void computeNeighbors();

    void interpolateLinear_CPU(int coarseProbeIdx, float3 pos, std::vector<float>& outDist, std::vector<float>& outDistSq);

    // ----------------------------------------------------------------
    // Resources & Debug
    // ----------------------------------------------------------------
    void uploadToGPU();
    void printDebugInfo(const std::string& filename);

    ref<Buffer> getProbeBuffer() const { return mpProbeBuffer; }
    ref<Buffer> getCornerBuffer() const { return mpCornerBuffer; }
    const std::vector<Probe>& getProbes() const { return mProbes; }
    // In your header (.h) or class definition:
    std::unordered_map<CornerKey, int, CornerKeyHasher> mCornerLookup;
    // ----------------------------------------------------------------
    // IO Interface
    // ----------------------------------------------------------------
    void saveToFile(const std::string& filename) const;
    void loadFromFile(const std::string& filename);

private:
    AdaptiveProbeVolume(ref<Device> pDevice);

    // ----------------------------------------------------------------
    // Internal Data Structures
    // ----------------------------------------------------------------

    ref<Device> mpDevice;

    std::vector<Probe> mProbes;
    std::vector<Corner> mCorners;

    // QUEUES
    std::vector<int> mPendingNewCorners;
    std::vector<int> mProbesPendingCheck;

    // Settings
    float mCurrentThreshold = 0.01f;
    int mMaxLevel = 7;
    bool mUseRelativeError = false;

    ref<Buffer> mpProbeBuffer;
    ref<Buffer> mpCornerBuffer;
};
