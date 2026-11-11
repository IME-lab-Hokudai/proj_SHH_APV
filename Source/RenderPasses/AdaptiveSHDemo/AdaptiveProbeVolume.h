#pragma once
#include <vector>
#include "AdaptiveGridData.slang"
#include "Falcor.h"

using namespace Falcor;

class AdaptiveProbeVolume : public Object
{
public:

    struct Probe
    {
        bool isLeaf = true;
        int level = 0;
        bool addedByEdgeAware = false;
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
        //std::vector<float> distMean;   // New
        //std::vector<float> distMeanSq; // New
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

    void interpolateHermite_CPU(int coarseProbeIdx, float3 pos, std::vector<float3>& outCoeffs, std::vector<GradSHCoeff>& outGrads);

    void constrainHangingNodesHermite();

    static ref<AdaptiveProbeVolume> create(ref<Device> pDevice);

    int traverseOctreeCPU(float3 pos) const;

    void computeNeighbors();

    // ----------------------------------------------------------------
    // Resources & Debug
    // ----------------------------------------------------------------
    void uploadToGPU();
    // Hanging-node constraint applied in uploadToGPU (seam option "constraint" of AdaptiveSHDemo; same rule as
    // AdaptiveProbeVolume::constrainHangingCorners in PrecomputeSHCoefficients): interp 0 Hermite, 1 trilinear,
    // 2 limited Hermite, -1 off. It changes the loaded corners, so load the grid again to switch it off.
    void setHangingConstraint(int interp) { mHangingConstraintInterp = interp; }
    uint32_t getConstrainedCornerCount() const { return mConstrainedCorners; }
    void printDebugInfo(const std::string& filename);

    ref<Buffer> getProbeBuffer() const { return mpProbeBuffer; }
    ref<Buffer> getCornerBuffer() const { return mpCornerBuffer; }
    const std::vector<Probe>& getProbes() const { return mProbes; }
    // In your header (.h) or class definition:
    std::unordered_map<CornerKey, int, CornerKeyHasher> mCornerLookup;

    void loadFromFile(const std::string& filename);
    ref<Buffer> getSeedProbeIndexBuffer() const { return mpSeedProbeIndexBuffer; }

    bool getUseSeedGrid() const { return mUseSeedGrid; }
    float3 getSeedMinPoint() const { return mSeedMinPoint; }
    float3 getSeedCellSize() const { return mSeedCellSize; }
    uint3 getSeedResolution() const { return mSeedResolution; }

    // Seed grid as uploaded (bind these to SeedGridCB): the file's seed grid, or for a grid without one a single
    // 1 x 1 x 1 cell covering the root, so the shader always starts from a seed cell and derives node bounds from it.
    float3 getGpuSeedMinPoint() const { return mGpuSeedMinPoint; }
    float3 getGpuSeedCellSize() const { return mGpuSeedCellSize; }
    uint3 getGpuSeedResolution() const { return mGpuSeedResolution; }
    // Size of the finest node (the lattice of the Oat seam handling in XAtlasSceneShader.slang, starting at the GPU seed
    // minimum); bind to SeedGridCB.gLatticeCellSize.
    float3 getGpuLatticeCellSize() const { return mGpuLatticeCellSize; }

    // Runtime GPU memory of the uploaded grid (compact layout), in bytes.
    uint64_t getGpuNodeBytes() const { return mGpuNodeBytes; }
    uint64_t getGpuCornerBytes() const { return mGpuCornerBytes; }
    uint64_t getGpuSeedBytes() const { return mGpuSeedBytes; }
private:
    AdaptiveProbeVolume(ref<Device> pDevice);

    ref<Device> mpDevice;

    std::vector<Probe> mProbes;
    std::vector<Corner> mCorners;

    // QUEUES
    std::vector<int> mPendingNewCorners;
    std::vector<int> mProbesPendingCheck;

    // Settings
    float mCurrentThreshold = 0.01f;
    int mMaxLevel = 5;
    int mHangingConstraintInterp = -1;
    uint32_t mConstrainedCorners = 0;
    bool mUseRelativeError = false;

    ref<Buffer> mpProbeBuffer;
    ref<Buffer> mpCornerBuffer;

    //statistic
    double mBuildTimeMs = 0.0;

    bool mUseSeedGrid = false;
    float3 mSeedMinPoint = float3(0.f);
    float3 mSeedCellSize = float3(0.f);
    uint3 mSeedResolution = uint3(0);
    std::vector<int> mSeedProbeIndices; // size = rx * ry * rz
    ref<Buffer> mpSeedProbeIndexBuffer;

    // Set by uploadToGPU.
    float3 mGpuSeedMinPoint = float3(0.f);
    float3 mGpuSeedCellSize = float3(0.f);
    uint3 mGpuSeedResolution = uint3(0);
    float3 mGpuLatticeCellSize = float3(1.f);
    uint64_t mGpuNodeBytes = 0, mGpuCornerBytes = 0, mGpuSeedBytes = 0;
};
