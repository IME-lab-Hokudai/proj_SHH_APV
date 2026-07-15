#pragma once
#include "Falcor.h"
#include <vector>
#include "envMap_SH.h"

using namespace Falcor;

class AdaptiveProbeVolume : public Object
{
public:
    struct PaperScalarStats
    {
        uint64_t count = 0;
        double sum = 0.0;
        double sumSq = 0.0;
        double minValue = std::numeric_limits<double>::infinity();
        double maxValue = -std::numeric_limits<double>::infinity();

        // Store values because median / p95 are useful in the paper.
        std::vector<float> values;

        void reset()
        {
            count = 0;
            sum = 0.0;
            sumSq = 0.0;
            minValue = std::numeric_limits<double>::infinity();
            maxValue = -std::numeric_limits<double>::infinity();
            values.clear();
        }

        void add(double v)
        {
            if (!std::isfinite(v))
                return;

            count++;
            sum += v;
            sumSq += v * v;
            minValue = std::min(minValue, v);
            maxValue = std::max(maxValue, v);
            values.push_back(float(v));
        }

        double mean() const
        {
            return count > 0 ? sum / double(count) : 0.0;
        }

        double variance() const
        {
            if (count <= 1)
                return 0.0;

            double m = mean();
            return std::max(0.0, sumSq / double(count) - m * m);
        }

        double stddev() const
        {
            return std::sqrt(variance());
        }

        double minOrZero() const
        {
            return count > 0 ? minValue : 0.0;
        }

        double maxOrZero() const
        {
            return count > 0 ? maxValue : 0.0;
        }
    };

    struct ResidualSamplePaperStats
    {
        uint64_t tracedCorners = 0;

        uint64_t noParentCorners = 0;
        uint64_t parentCorners = 0;

        uint64_t validResidualCorners = 0;
        uint64_t invalidResidualCorners = 0;

        uint64_t scaleBelowOne = 0;
        uint64_t scaleEqualOne = 0;
        uint64_t scaleAboveOne = 0;

        uint64_t scaleAtMinClamp = 0;
        uint64_t scaleAtMaxClamp = 0;

        PaperScalarStats observedResidual;
        PaperScalarStats predictedParentError;
        PaperScalarStats residualRatio;
        PaperScalarStats correctionScale;

        void reset()
        {
            tracedCorners = 0;

            noParentCorners = 0;
            parentCorners = 0;

            validResidualCorners = 0;
            invalidResidualCorners = 0;

            scaleBelowOne = 0;
            scaleEqualOne = 0;
            scaleAboveOne = 0;

            scaleAtMinClamp = 0;
            scaleAtMaxClamp = 0;

            observedResidual.reset();
            predictedParentError.reset();
            residualRatio.reset();
            correctionScale.reset();
        }
    };

    struct ResidualDecisionPaperStats
    {
        uint64_t testedCells = 0;

        // Both methods make same decision.
        uint64_t bothSubdivide = 0;
        uint64_t bothLeaf = 0;

        // Important paper statistic:
        // original Hessian would subdivide,
        // residual-calibrated metric keeps the cell as leaf.
        uint64_t prunedByResidual = 0;

        // Opposite case:
        // original Hessian would keep leaf,
        // residual-calibrated metric subdivides.
        uint64_t addedByResidual = 0;

        PaperScalarStats originalHessianError;
        PaperScalarStats correctedResidualError;
        PaperScalarStats correctedOverOriginalRatio;

        void reset()
        {
            testedCells = 0;

            bothSubdivide = 0;
            bothLeaf = 0;

            prunedByResidual = 0;
            addedByResidual = 0;

            originalHessianError.reset();
            correctedResidualError.reset();
            correctedOverOriginalRatio.reset();
        }
    };

    static constexpr int kPaperStatsMaxLevel = 16;

    struct ResidualPaperStats
    {
        ResidualSamplePaperStats allSamples;
        ResidualDecisionPaperStats allDecisions;

        std::array<ResidualSamplePaperStats, kPaperStatsMaxLevel> samplesByParentLevel;
        std::array<ResidualDecisionPaperStats, kPaperStatsMaxLevel> decisionsByCellLevel;

        void reset()
        {
            allSamples.reset();
            allDecisions.reset();

            for (auto& s : samplesByParentLevel)
                s.reset();

            for (auto& d : decisionsByCellLevel)
                d.reset();
        }
    };

    void resetResidualPaperStats();

    void printResidualPaperStats(std::ostream& out) const;

    void exportResidualPaperStatsCSV(
        const std::string& summaryCsv,
        const std::string& levelCsv
    ) const;

    void writeResidualPaperStatsLog(const std::string& filename) const;

    void printResidualMainStats(std::ostream& out) const;

    void writeResidualMainStatsLog(
        const std::string& filename
    ) const;

    void exportResidualMainStatsCSV(
        const std::string& filename
    ) const;


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

        float lastError = 0.0f;
        bool needsRecheck = false;
    };
    struct Corner
    {
        float3 position;

        // Physics Data (Full bands)
        std::vector<float3> shCoeffs;
        std::vector<GradSHCoeff> shGradients;

        float maxLambdaVecL2 = 0.0f; // Curvature
        float coeffVecL2 = 0.0f; // ||L||

         // ------------------------------------------------------------
        // Residual-corrected Hessian scale
        // ------------------------------------------------------------
        // Parent cell that created this corner.
        // -1 means this corner has no residual correction.
        int residualParentProbeIdx = -1;

        // Hessian-predicted parent error used as denominator.
        // This corresponds to parent-level E_abs prediction.
        float residualPredictedError = 0.0f;

        // Observed parent interpolation residual:
        // R_abs = ||f_actual(x) - f_parent_interp(x)||
        float residualObservedError = 0.0f;

        // rho = R_abs / (predicted E_abs + eps)
        float residualRatio = 1.0f;

        // Final correction scale:
        // s = clamp((1 - eta) + eta * rho, minScale, maxScale)
        float residualCorrectionScale = 1.0f;
    };

    struct MemoryFootprintInfo
    {
        uint64_t gpuCornersBytes = 0;
        uint64_t gpuProbesBytes = 0;
        uint64_t totalBytes = 0;
        float totalMB = 0.0; 
    };

    struct BuildStats
    {
        uint64_t totalEvaluatedCorners = 0;
        uint64_t totalRayCount = 0;
        uint64_t dispatchCount = 0;
        uint64_t readbackCount = 0;

        uint64_t buildGpuBytes = 0;
        double buildGpuMB = 0.0;

        uint64_t runtimeGpuBytes = 0;
        double runtimeGpuMB = 0.0;
    };

    struct ProgressiveStage
    {
        uint32_t spp;
        float thresholdScale;
        bool storeFullRuntimeData;
        bool refineAfterThisStage;
    };

    std::vector<ProgressiveStage> stages =
    {
        {128,  0.50f, false, true},
        {512,  0.80f, false, true},
        {1024, 1.00f, true,  false}
    };

    MemoryFootprintInfo calculateMemoryFootprint() const;
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

    enum class ErrorMetricMode
    {
        SHSpace,
        IrradianceSpace
    };

    void setErrorMetricMode(ErrorMetricMode mode) { mErrorMetricMode = mode; }

    void setResidualCorrection(
        bool enabled,
        float eta = 1.0f,
        float minScale = 0.5f,
        float maxScale = 2.0f
    )
    {
        mUseResidualCorrection = enabled;
        mResidualCorrectionEta = eta;
        mResidualCorrectionMinScale = minScale;
        mResidualCorrectionMaxScale = maxScale;
    }

    static ref<AdaptiveProbeVolume> create(ref<Device> pDevice);

    void startBuild(const ref<Scene>& pScene, float errorThreshold, bool useRelativeError = false);
    void interpolateHermite_CPU(int coarseProbeIdx, float3 pos, std::vector<float3>& outCoeffs, std::vector<GradSHCoeff>& outGrads);
    void interpolateHermite_CPU(int coarseProbeIdx, float3 pos, std::vector<float3>& outCoeffs) const;
    int getProbeLevelCPU(int probeIdx) const
    {
        if (probeIdx < 0 || probeIdx >= int(mProbes.size()))
            return -1;

        return mProbes[probeIdx].level;
    }
    void constrainHangingNodesHermite();
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
        const std::vector<float3x3>& hessians
        //const std::vector<float>& distMeans,
        //const std::vector<float>& distMeanSqs
    );

    // Calculate Error per Probe, Subdivide if needed, Generate new Corners
    void finishBatch();

    int traverseOctreeCPU(float3 pos) const;

    void computeNeighbors();

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
    double getBuildTimeMs() const { return mBuildTimeMs; }
    void setBuildTimeMs(double t) { mBuildTimeMs = t; }

    void startBuildSeeded(
        const ref<Scene>& pScene,
        uint3 seedResolution,
        float errorThreshold,
        bool useRelativeError = false
    );

    void getPendingPositionsRange(
        uint32_t start,
        uint32_t count,
        std::vector<float3>& positions
    ) const;

    void setCornerDataRange(
        uint32_t start,
        const std::vector<std::vector<float3>>& coeffsBatch,
        const std::vector<std::vector<GradSHCoeff>>& gradsBatch,
        const std::vector<std::vector<float3x3>>& hessiansBatch
    );

    float3 evaluateIrradianceHermiteCPU(float3 posW, float3 normalW) const;

    uint32_t getPendingCornerCount() const { return (uint32_t)mPendingNewCorners.size(); }

    // New: seed-grid aware traversal helpers
    int findSeedProbeCPU(float3 pos) const;
    bool hasSeedGrid() const { return mUseSeedGrid; }

    // Optional helper if you seed externally
    void setSeedGridMetadata(float3 minPoint, float3 cellSize, uint3 resolution, const std::vector<int>& seedProbeIndices)
    {
        mUseSeedGrid = true;
        mSeedMinPoint = minPoint;
        mSeedCellSize = cellSize;
        mSeedResolution = resolution;
        mSeedProbeIndices = seedProbeIndices;
    }

    //==== progressive build ====//

    const BuildStats& getBuildStats() const { return mBuildStats; }
    void resetBuildStats()
    {
        mBuildStats = BuildStats{};
    }

    void recordTraceBatch(uint32_t cornerCount, uint32_t spp)
    {
        mBuildStats.totalEvaluatedCorners += cornerCount;
        mBuildStats.totalRayCount += uint64_t(cornerCount) * uint64_t(spp);
        mBuildStats.dispatchCount++;
        mBuildStats.readbackCount++;
    }

    void setCornerMetricData(
        uint32_t batchIndex,
        float coeffVecL2,
        float maxLambdaVecL2
    );

    void scheduleLeafCornersForRefinementRecheck();
    void setCurrentThreshold(float threshold) { mCurrentThreshold = threshold; };
    void setCornerRuntimeData(
        uint32_t batchIndex,
        const std::vector<float3>& coeffs,
        const std::vector<GradSHCoeff>& grads
    );

    void clearPendingBatch();
    void scheduleAllLeafCornersForRuntimeBake();
    void printCoarseStageDebugInfo(const std::string& filename) const;
    void finishBatchCoarseLimited(int maxEvalLevel);
private:
    AdaptiveProbeVolume(ref<Device> pDevice);
    BuildStats mBuildStats;
    // ----------------------------------------------------------------
    // Internal Data Structures
    // ----------------------------------------------------------------
    ErrorMetricMode mErrorMetricMode = ErrorMetricMode::SHSpace;
    ref<Device> mpDevice;

    std::vector<Probe> mProbes;
    std::vector<Corner> mCorners;

    // QUEUES
    std::vector<int> mPendingNewCorners;
    std::vector<int> mProbesPendingCheck;

    // Settings
    float mCurrentThreshold = 0.01f;
    //int mMaxLevel = 6;
    int mMaxLevel = 5;
    //int mMaxLevel = 0;
    bool mUseRelativeError = false;

    ref<Buffer> mpProbeBuffer;
    ref<Buffer> mpCornerBuffer;

    //statistic
    double mBuildTimeMs = 0.0;

    // Top-level coarse seed grid metadata
    bool mUseSeedGrid = false;
    float3 mSeedMinPoint = float3(0.f);
    float3 mSeedCellSize = float3(0.f);
    uint3 mSeedResolution = uint3(0);
    std::vector<int> mSeedProbeIndices; // size = rx * ry * rz

    bool mUseResidualCorrection = true;
    float mResidualCorrectionEta = 1.0f;
    float mResidualCorrectionMinScale = 0.5f;
    float mResidualCorrectionMaxScale = 2.0f;
    float mResidualCorrectionEps = 1e-6f;
    float computeParentPointHessianPrediction(
        int parentProbeIdx,
        const float3& position
    ) const;

private:
        ResidualPaperStats mResidualPaperStats;

        void recordResidualSampleForPaper(
            int parentLevel,
            bool hasParent,
            bool validResidual,
            float observedResidual,
            float predictedParentError,
            float rho,
            float scale
        );

        void recordResidualDecisionForPaper(
            int cellLevel,
            float originalHessianError,
            float correctedResidualError
        );
};
