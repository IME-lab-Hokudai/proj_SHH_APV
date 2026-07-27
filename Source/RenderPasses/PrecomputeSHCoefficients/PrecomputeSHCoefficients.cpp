/***************************************************************************
 # Copyright (c) 2015-23, NVIDIA CORPORATION. All rights reserved.
 #
 # Redistribution and use in source and binary forms, with or without
 # modification, are permitted provided that the following conditions
 # are met:
 #  * Redistributions of source code must retain the above copyright
 #    notice, this list of conditions and the following disclaimer.
 #  * Redistributions in binary form must reproduce the above copyright
 #    notice, this list of conditions and the following disclaimer in the
 #    documentation and/or other materials provided with the distribution.
 #  * Neither the name of NVIDIA CORPORATION nor the names of its
 #    contributors may be used to endorse or promote products derived
 #    from this software without specific prior written permission.
 #
 # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS "AS IS" AND ANY
 # EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 # IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 # PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 # CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 # EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 # PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 # PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 # OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 # (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 # OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 **************************************************************************/
#define PROBE_MODE_ADAPTIVE 0
#define PROBE_MODE_UNIFORM  1
 // CHANGE THIS LINE TO SWITCH MODES:
//#define CURRENT_PROBE_MODE PROBE_MODE_UNIFORM
#define CURRENT_PROBE_MODE PROBE_MODE_ADAPTIVE

#include <fstream>
#include "PrecomputeSHCoefficients.h"

#include "envMap_SH.h"
#include "Rendering/Lights/EmissivePowerSampler.h"
#include "Rendering/Lights/EmissiveUniformSampler.h"
#include "Rendering/Lights/LightBVHSampler.h"
#include <cmath>
#include "ProbeSamplingData.slang"
#include <Scene/Material/StandardMaterial.h>
#include <chrono>
const int numSamplesPerProbe = 4096;
//const int numSamplesPerProbe = 1024;
//const int numSamplesPerProbe = 2048;
const uint32_t kMaxSamplesPerProbe = 1024; //used in abandoned progressive build test.

//const int numSamplesPerProbe = 1;
const int verificationRes = 100;
const float verificationH = 0.001f;
const float verificationY = 0.2f;
const float verificationExtent = 0.25f;

//const float ErrorThreshold = 100000000.0f;
//const float ErrorThreshold = 20.0f;
//const float ErrorThreshold = 10.0f;
//const float ErrorThreshold = 8.5f;
//const float ErrorThreshold = 8.0f;
//const float ErrorThreshold = 7.95f;
//const float ErrorThreshold = 4.25f;
//const float ErrorThreshold = 3.0f;
//const float ErrorThreshold = 5.0f;
//const float ErrorThreshold = 1.0f;
//const float ErrorThreshold =3.5f;//threshold for Erel
//const float ErrorThreshold =1.5f;//threshold for Erel
const float ErrorThreshold =2.0f;//threshold for Erel
const bool useRelativeError = false;
const bool useIrradianceSpaceMetric = false;
const bool useResidualCorrection = true;
//const bool useResidualCorrection = false;

const float residualPruneStrength = 0.00f;
const float residualRefineStrength = 0.50f;

const float residualCorrectionMinScale = 0.5f;
const float residualCorrectionMaxScale = 2.0f;

const float residualConfidenceTau0 = 0.9f;
const float residualConfidenceTau1 = 1.1f;
const float residualConfidenceEps = 1e-3f;
//const bool useIrradianceSpaceMetric = false;
//const bool useRelativeError = true;
//const uint3 unifromGridSize = uint3(16, 16, 16);
const uint3 unifromGridSize = uint3(32, 32, 32);
//const uint3 unifromGridSize = uint3(8, 8, 8);
//const uint3 unifromGridSize = uint3(64, 64, 64);
//const std::string loadFromFileName = "DirectAbsErr8p5N6DataScene.txt";
const std::string loadFromFileName = "DirectAbsErr2HessianMetricCornellThinSlabV2.txt";

//const std::string saveToFileName = "Seeded16DirectAbsErr3SubwayCorridorNoOpen.txt";
//const std::string saveToFileName = "Seeded16DirectAbsErr5SubwayCorridorNoOpen.txt";
//const std::string saveToFileName = "Seeded16DirectAbsErr10SubwayCorridorNoOpen.txt";

//const std::string saveToFileName = "Seeded8DirectAbsErr5SubwayCorridorNoOpen.txt";
//const std::string saveToFileName = "Seeded8RelErr3SubwayCorridorNoOpen.txt";
//const std::string saveToFileName = "Seeded8RelErr2SubwayCorridorNoOpen.txt";
//const std::string saveToFileName = "Seeded8NewMetricAbsErr5SubwayCorridorNoOpen.txt";

//const std::string saveToFileName = "DirectUniformGrid32SubwayCorridorNoOpen.txt";
//const std::string saveToFileName = "DirectUniformGrid64SubwayCorridorNoOpen.txt";

//const std::string saveToFileName = "IndirectUniformGrid64SubwayCorridorNoOpen.txt";
//const std::string saveToFileName = "IndirectUniformGrid32SubwayCorridorNoOpen.txt";

//const std::string saveToFileName = "DirectTestDC.txt";
//const std::string saveToFileName = "DirectAbsErr10AsymScene.txt";
//const std::string saveToFileName = "DirectAbsErr10AsymSceneN5.txt";
//const std::string saveToFileName = "DirectAbsErr10AsymSceneN6.txt";
//const std::string saveToFileName = "DirectAbsErr10AsymSceneN6Progressive.txt";
//const std::string saveToFileName = "DirectU32AsymScene.txt";

//const std::string saveToFileName = "DirectAbsErr20DataScene.txt";
//const std::string saveToFileName = "DirectAbsErr10DataScene.txt";
//const std::string saveToFileName = "DirectAbsErr5DataScene.txt";
//const std::string saveToFileName = "DirectAbsErr1DataScene.txt";

//const std::string saveToFileName = "U64DataScene.txt";
//const std::string saveToFileName = "U32DataScene.txt";
//const std::string saveToFileName = "U32DataScene_4096spp_4spt.txt";

//const std::string saveToFileName = "DirectAbsErr20N6DataScene.txt";
//const std::string saveToFileName = "DirectAbsErr10N6DataScene.txt";
//const std::string saveToFileName = "DirectAbsErr8p5N6DataScene.txt";
//const std::string saveToFileName = "DirectAbsErr8N6DataScene.txt";
//const std::string saveToFileName = "DirectAbsErr5DataScene.txt";
//const std::string saveToFileName = "DirectAbsErr1DataScene.txt";

//const std::string saveToFileName = "DirectAbsErr8p5N6IndirectDataScene.txt";

//const std::string saveToFileName = "DirectAbsErr10DataScenePresentativeNormalMetric.txt";
//const std::string saveToFileName = "DirectAbsErr8DataSceneIrrMetric.txt";
//const std::string saveToFileName = "DirectAbsErr8DataScenePresentativeNormalMetricAvg.txt";
//const std::string saveToFileName = "DirectAbsErr3DataScenePresentativeNormalMetricAvg.txt";
//const std::string saveToFileName = "DirectAbsErr2p1DataScenePresentativeNormalMetricAvg.txt";
//const std::string saveToFileName = "DirectAbsErr2DataScenePresentativeNormalMetricAvg.txt";
//const std::string saveToFileName = "DirectAbsErr8DataScenePresentativeNormalMetric.txt";
//const std::string saveToFileName = "DirectAbsErr7p5DataSceneIrrMetric.txt";
//const std::string saveToFileName = "DirectAbsErr5DataSceneIrrMetric.txt";
//const std::string saveToFileName = "DirectAbsErr4p25DataSceneIrrMetric.txt";
//const std::string saveToFileName = "DirectAbsErr4DataSceneIrrMetric.txt";
//const std::string saveToFileName = "DirectAbsErr3DataSceneIrrMetric.txt";
//const std::string saveToFileName = "DirectAbsErr2p1DataSceneIrrMetric.txt";
//const std::string saveToFileName = "DirectAbsErr20DataSceneAvg.txt";
//const std::string saveToFileName = "DirectAbsErr10DataSceneAvg.txt";
//const std::string saveToFileName = "DirectAbsErr5DataSceneAvg.txt";
//const std::string saveToFileName = "DirectAbsErr4DataSceneAvg.txt";
//const std::string saveToFileName = "DirectAbsErr3p5DataSceneAvg.txt";
//const std::string saveToFileName = "DirectAbsErr3DataSceneAvg.txt";
//const std::string saveToFileName = "DirectAbsErr2DataSceneAvg.txt";

//const std::string saveToFileName = "DirectAbsErr10ResidualScaleMetric.txt";
//const std::string saveToFileName = "DirectAbsErr10ResidualScaleV2Metric.txt";
//const std::string saveToFileName = "DirectAbsErr5ResidualScaleV2Metric.txt";
//const std::string saveToFileName = "DirectAbsErr5ResidualScaleV2MetricCornellThinSlab.txt";
//const std::string saveToFileName = "DirectAbsErr2CornellThinSlab.txt";
//const std::string saveToFileName = "DirectAbsErr2ResidualScaleV2MetricCornellThinSlab.txt";
//const std::string saveToFileName = "DirectAbsErr5ResidualScaleMetric.txt";
//const std::string saveToFileName = "DirectAbsErr2HessianMetricCornellThinSlab.txt";
//const std::string saveToFileName = "DirectAbsErr2HessianMetricCornellThinSlabV2.txt";
const std::string saveToFileName = "DirectAbsErr2EdgeMetricCornellThinSlabV2.txt";

//const std::string saveToFileName = "U64CornellShadowBoundaryScene.txt";
//const std::string saveToFileName = "U32CornellShadowBoundaryScene.txt";
namespace
{
//const char kShaderFile[] = "RenderPasses/PrecomputeSHCoefficients/SHShader.slang";
#if CURRENT_PROBE_MODE == PROBE_MODE_ADAPTIVE
    const char kShaderFile[] = "RenderPasses/PrecomputeSHCoefficients/SHAdaptiveProbeShader.slang";
#else
    const char kShaderFile[] = "RenderPasses/PrecomputeSHCoefficients/SHGridShader.slang";
#endif
const char kEnvMapShaderFile[] = "RenderPasses/PrecomputeSHCoefficients/EnvMapShader.slang";
const char kProbeSamplingFile[] = "RenderPasses/PrecomputeSHCoefficients/ProbeSampling.rt.slang";
const char kShowReconstructedEnvMap[] = "Show environment map";
const char kShowSHGrid[] = "Show SH grid";
const char kShowEdgeAddedVoxels[] ="Show edge-added voxels";
} // namespace

namespace
{
    static std::string replaceExtensionOrAppend(
        const std::string& baseName,
        const std::string& suffixWithExtension
    )
    {
        std::string result = baseName;

        size_t dotPos = result.find_last_of('.');

        if (dotPos != std::string::npos)
        {
            result.replace(
                dotPos,
                std::string::npos,
                suffixWithExtension
            );
        }
        else
        {
            result += suffixWithExtension;
        }

        return result;
    }

    static float luminance709(float3 c)
    {
        return dot(c, float3(0.2126f, 0.7152f, 0.0722f));
    }
    static double percentileVector(std::vector<float> values, double percentile)
    {
        if (values.empty())
            return 0.0;

        std::sort(values.begin(), values.end());

        double pos = percentile * 0.01 * double(values.size() - 1);
        size_t i0 = size_t(std::floor(pos));
        size_t i1 = std::min(i0 + 1, values.size() - 1);

        double t = pos - double(i0);

        return double(values[i0]) * (1.0 - t) + double(values[i1]) * t;
    }

    struct DeltaErrorStats
    {
        std::vector<float> deltaAPEValues;

        double sumDeltaAPE = 0.0;
        double sumDeltaAbsErr = 0.0;

        uint64_t sampleCount = 0;
        uint64_t adaptiveBetterCount = 0;
    };

    static void addDeltaErrorSample(
        DeltaErrorStats& stats,
        float3 u32Irr,
        float3 adaptiveIrr,
        float3 refIrr,
        float relEpsilon
    )
    {
        float yRef = luminance709(refIrr);
        float yU32 = luminance709(u32Irr);
        float yAdaptive = luminance709(adaptiveIrr);

        float absErrU32 = std::abs(yU32 - yRef);
        float absErrAdaptive = std::abs(yAdaptive - yRef);

        float denom = std::max(std::abs(yRef), relEpsilon);

        float apeU32 = absErrU32 / denom * 100.0f;
        float apeAdaptive = absErrAdaptive / denom * 100.0f;

        float deltaAPE = apeU32 - apeAdaptive;
        float deltaAbsErr = absErrU32 - absErrAdaptive;

        if (!std::isfinite(deltaAPE))
            return;

        stats.deltaAPEValues.push_back(deltaAPE);
        stats.sumDeltaAPE += double(deltaAPE);
        stats.sumDeltaAbsErr += double(deltaAbsErr);
        stats.sampleCount++;

        if (absErrAdaptive < absErrU32)
            stats.adaptiveBetterCount++;
    }

    static void writeDeltaErrorRow(
        std::ofstream& csv,
        const std::string& region,
        const DeltaErrorStats& stats,
        uint64_t totalSampleCount
    )
    {
        double sampleSharePercent = totalSampleCount > 0
            ? 100.0 * double(stats.sampleCount) / double(totalSampleCount)
            : 0.0;

        double meanDeltaMAPE = stats.sampleCount > 0
            ? stats.sumDeltaAPE / double(stats.sampleCount)
            : 0.0;

        double meanDeltaAbsErr = stats.sampleCount > 0
            ? stats.sumDeltaAbsErr / double(stats.sampleCount)
            : 0.0;

        double medianDeltaMAPE = percentileVector(stats.deltaAPEValues, 50.0);
        double p05DeltaMAPE = percentileVector(stats.deltaAPEValues, 5.0);
        double p95DeltaMAPE = percentileVector(stats.deltaAPEValues, 95.0);

        double adaptiveBetterRate = stats.sampleCount > 0
            ? 100.0 * double(stats.adaptiveBetterCount) / double(stats.sampleCount)
            : 0.0;

        csv
            << region << ","
            << stats.sampleCount << ","
            << sampleSharePercent << ","
            << meanDeltaMAPE << ","
            << medianDeltaMAPE << ","
            << p05DeltaMAPE << ","
            << p95DeltaMAPE << ","
            << meanDeltaAbsErr << ","
            << adaptiveBetterRate
            << "\n";
    }

    struct IrradianceFieldErrorStats
    {
        std::vector<float> apeValues;

        double sumAbsErr = 0.0;
        double sumSqErr = 0.0;
        double sumRefLum = 0.0;
        double sumTestLum = 0.0;

        uint64_t sampleCount = 0;
    };

    static void addIrradianceErrorSample(
        IrradianceFieldErrorStats& stats,
        float3 testIrr,
        float3 refIrr,
        float relEpsilon
    )
    {
        float yTest = luminance709(testIrr);
        float yRef = luminance709(refIrr);

        float absErr = std::abs(yTest - yRef);

        // Use abs(yRef) because SH irradiance can become slightly negative.
        float denom = std::max(std::abs(yRef), relEpsilon);
        float ape = absErr / denom * 100.0f;

        if (!std::isfinite(ape))
            return;

        stats.apeValues.push_back(ape);
        stats.sumAbsErr += double(absErr);
        stats.sumSqErr += double(absErr) * double(absErr);
        stats.sumRefLum += double(yRef);
        stats.sumTestLum += double(yTest);
        stats.sampleCount++;
    }

    static double meanVector(const std::vector<float>& values)
    {
        if (values.empty())
            return 0.0;

        double sum = 0.0;
        for (float v : values)
            sum += double(v);

        return sum / double(values.size());
    }



    static std::vector<float3> makeSixAxisNormalsWorld()
    {
        // IMPORTANT:
        // These are Falcor/world-space normals.
        // Do NOT pre-swizzle them. The CPU irradiance evaluator should mirror
        // the shader and do x=n.z, y=n.x, z=n.y internally.
        return {
            float3(1.0f,  0.0f,  0.0f),
            float3(-1.0f,  0.0f,  0.0f),
            float3(0.0f,  1.0f,  0.0f),
            float3(0.0f, -1.0f,  0.0f),
            float3(0.0f,  0.0f,  1.0f),
            float3(0.0f,  0.0f, -1.0f),
        };
    }
    struct LocalMethodErrorStats
    {
        double sumAPE = 0.0;
        double sumAbsError = 0.0;
        double sumSquaredError = 0.0;

        uint64_t sampleCount = 0;

        void add(float testLuminance, float referenceLuminance, float relEpsilon)
        {
            const float error =
                testLuminance - referenceLuminance;

            const float absError =
                std::abs(error);

            const float denominator =
                std::max(std::abs(referenceLuminance), relEpsilon);

            const float ape =
                absError / denominator * 100.0f;

            if (!std::isfinite(ape) ||
                !std::isfinite(absError))
            {
                return;
            }

            sumAPE += double(ape);
            sumAbsError += double(absError);
            sumSquaredError +=
                double(error) * double(error);

            sampleCount++;
        }

        double getMAPE() const
        {
            return sampleCount > 0
                ? sumAPE / double(sampleCount)
                : 0.0;
        }

        double getMAE() const
        {
            return sampleCount > 0
                ? sumAbsError / double(sampleCount)
                : 0.0;
        }

        double getRMSE() const
        {
            return sampleCount > 0
                ? std::sqrt(
                    sumSquaredError /
                    double(sampleCount)
                )
                : 0.0;
        }

        double getPSNR(float peakValue = 1.0f) const
        {
            const double rmse = getRMSE();

            if (rmse <= 0.0)
                return std::numeric_limits<double>::infinity();

            return 20.0 * std::log10(
                double(peakValue) / rmse
            );
        }
    };


    struct LocalEGCComparisonStats
    {
        LocalMethodErrorStats hessian;
        LocalMethodErrorStats egc;

        uint64_t egcBetterCount = 0;
        uint64_t hessianBetterCount = 0;
        uint64_t equalCount = 0;

        void add(
            float3 hessianIrradiance,
            float3 egcIrradiance,
            float3 referenceIrradiance,
            float relEpsilon
        )
        {
            const float referenceLuminance =
                luminance709(referenceIrradiance);

            const float hessianLuminance =
                luminance709(hessianIrradiance);

            const float egcLuminance =
                luminance709(egcIrradiance);

            const float hessianAbsError =
                std::abs(
                    hessianLuminance -
                    referenceLuminance
                );

            const float egcAbsError =
                std::abs(
                    egcLuminance -
                    referenceLuminance
                );

            hessian.add(
                hessianLuminance,
                referenceLuminance,
                relEpsilon
            );

            egc.add(
                egcLuminance,
                referenceLuminance,
                relEpsilon
            );

            constexpr float equalEpsilon = 1e-8f;

            if (egcAbsError + equalEpsilon < hessianAbsError)
            {
                egcBetterCount++;
            }
            else if (
                hessianAbsError + equalEpsilon <
                egcAbsError
                )
            {
                hessianBetterCount++;
            }
            else
            {
                equalCount++;
            }
        }

        uint64_t getComparisonCount() const
        {
            return
                egcBetterCount +
                hessianBetterCount +
                equalCount;
        }

        double getEGCBetterRate() const
        {
            const uint64_t count =
                getComparisonCount();

            return count > 0
                ? 100.0 *
                double(egcBetterCount) /
                double(count)
                : 0.0;
        }
    };

    struct SampleGainLossStats
    {
        // EGC produces lower absolute error than Hessian.
        uint64_t improvedCount = 0;

        // EGC produces higher absolute error than Hessian.
        uint64_t worsenedCount = 0;

        // Difference is within the numerical tolerance.
        uint64_t comparableCount = 0;

        uint64_t invalidCount = 0;

        // Positive magnitudes only.
        double totalGain = 0.0;
        double totalLoss = 0.0;

        std::vector<float> gainValues;
        std::vector<float> lossValues;

        void add(
            float hessianAbsError,
            float egcAbsError
        )
        {
            constexpr float kComparisonEps = 1e-6f;

            if (!std::isfinite(hessianAbsError) ||
                !std::isfinite(egcAbsError))
            {
                invalidCount++;
                return;
            }

            // Positive means EGC removed error.
            // Negative means EGC introduced error.
            const float errorReduction =
                hessianAbsError - egcAbsError;

            if (errorReduction > kComparisonEps)
            {
                improvedCount++;

                totalGain += double(errorReduction);

                gainValues.push_back(
                    errorReduction
                );
            }
            else if (errorReduction < -kComparisonEps)
            {
                const float loss =
                    -errorReduction;

                worsenedCount++;

                totalLoss += double(loss);

                lossValues.push_back(
                    loss
                );
            }
            else
            {
                comparableCount++;
            }
        }

        uint64_t getValidSampleCount() const
        {
            return
                improvedCount +
                worsenedCount +
                comparableCount;
        }

        double percent(uint64_t count) const
        {
            const uint64_t validCount =
                getValidSampleCount();

            return validCount > 0
                ? 100.0 *
                double(count) /
                double(validCount)
                : 0.0;
        }

        double getMeanGain() const
        {
            return improvedCount > 0
                ? totalGain /
                double(improvedCount)
                : 0.0;
        }

        double getMedianGain() const
        {
            return percentileVector(
                gainValues,
                50.0
            );
        }

        double getMeanLoss() const
        {
            return worsenedCount > 0
                ? totalLoss /
                double(worsenedCount)
                : 0.0;
        }

        double getMedianLoss() const
        {
            return percentileVector(
                lossValues,
                50.0
            );
        }

        double getNetGain() const
        {
            return totalGain - totalLoss;
        }

        double getNetMeanErrorReduction() const
        {
            const uint64_t validCount =
                getValidSampleCount();

            return validCount > 0
                ? getNetGain() /
                double(validCount)
                : 0.0;
        }

        double getGainLossRatio() const
        {
            if (totalLoss > 0.0)
                return totalGain / totalLoss;

            if (totalGain > 0.0)
                return std::numeric_limits<double>::infinity();

            return 0.0;
        }
    };

    static void writeLocalEGCComparisonRow(
        std::ofstream& csv,
        const std::string& regionName,
        int probeIndex,
        int probeLevel,
        float3 minPoint,
        float3 maxPoint,
        const LocalEGCComparisonStats& stats
    )
    {
        csv << std::fixed << std::setprecision(8);

        csv
            << regionName << ","
            << probeIndex << ","
            << probeLevel << ","

            << minPoint.x << ","
            << minPoint.y << ","
            << minPoint.z << ","

            << maxPoint.x << ","
            << maxPoint.y << ","
            << maxPoint.z << ","

            << stats.hessian.sampleCount << ","

            << stats.hessian.getMAPE() << ","
            << stats.egc.getMAPE() << ","
            << stats.hessian.getMAPE() -
            stats.egc.getMAPE() << ","

            << stats.hessian.getMAE() << ","
            << stats.egc.getMAE() << ","
            << stats.hessian.getMAE() -
            stats.egc.getMAE() << ","

            << stats.hessian.getRMSE() << ","
            << stats.egc.getRMSE() << ","

            << stats.hessian.getPSNR() << ","
            << stats.egc.getPSNR() << ","

            << stats.egcBetterCount << ","
            << stats.hessianBetterCount << ","
            << stats.equalCount << ","
            << stats.getEGCBetterRate()
            << "\n";
    }
    static void writeFieldErrorSummaryRow(
        std::ofstream& csv,
        const std::string& method,
        const std::string& fileName,
        int maxDepth,
        const std::string& threshold,
        const IrradianceFieldErrorStats& stats
    )
    {
        double mape = meanVector(stats.apeValues);
        double p50 = percentileVector(stats.apeValues, 50.0);
        double p95 = percentileVector(stats.apeValues, 95.0);
        double p99 = percentileVector(stats.apeValues, 99.0);

        double mae = stats.sampleCount > 0
            ? stats.sumAbsErr / double(stats.sampleCount)
            : 0.0;

        double rmse = stats.sampleCount > 0
            ? std::sqrt(stats.sumSqErr / double(stats.sampleCount))
            : 0.0;

        double meanRef = stats.sampleCount > 0
            ? stats.sumRefLum / double(stats.sampleCount)
            : 0.0;

        double meanTest = stats.sampleCount > 0
            ? stats.sumTestLum / double(stats.sampleCount)
            : 0.0;

        csv
            << method << ","
            << fileName << ","
            << maxDepth << ","
            << threshold << ","
            << stats.sampleCount << ","
            << mape << ","
            << p50 << ","
            << p95 << ","
            << p99 << ","
            << mae << ","
            << rmse << ","
            << meanRef << ","
            << meanTest
            << "\n";
    }

    static void writeFieldErrorSamples(
        std::ofstream& csv,
        const std::string& method,
        const IrradianceFieldErrorStats& stats,
        uint32_t stride
    )
    {
        stride = std::max(1u, stride);

        for (size_t i = 0; i < stats.apeValues.size(); i += stride)
        {
            csv << method << "," << stats.apeValues[i] << "\n";
        }
    }

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

    //error vs distance
    static float maxAbsEigenvalue3x3(const float3x3& H)
    {
        float e1, e2, e3;
        computeEigenvalues3x3(H, e1, e2, e3);
        return std::max({ std::abs(e1), std::abs(e2), std::abs(e3) });
    }

    static float lambdaL2OverBands(const std::array<float3x3, 9>& H)
    {
        float sumSq = 0.0f;

        for (int basisIdx = 0; basisIdx < 9; ++basisIdx)
        {
            float lambda = maxAbsEigenvalue3x3(H[basisIdx]);
            sumSq += lambda * lambda;
        }

        return std::sqrt(sumSq);
    }

    void exportAnalyticErrorVsDistanceTest()
    {
        const std::string outCsv = "AnalyticErrorVsDistanceTest.csv";

        // Virtual grid only defines ||Delta x|| for Eabs.
        const float sceneExtent = 1.0f;
        const uint32_t virtualGridResolution = 64;
        const float cellSize = sceneExtent / float(virtualGridResolution);
        const float deltaXNorm = std::sqrt(3.0f) * cellSize;

        // Controlled setup:
        // fixed surface point s = origin.
        // q = s - X = r * qDir.
        const float3 qDir = math::normalize(float3(0.35f, 0.45f, 0.82f));

        // Surface normal faces the probe direction.
        const float3 n = -qDir;

        // Constant radiance assumption.
        const float L = 1.0f;

        // One analytic patch contribution.
        // Same role as Omega_i = 4pi / sampleCount in your normal code.
        const float Omega_i = 1.0f;

        const float rMin = 0.02f;
        const float rMax = 4.0f;
        const uint32_t sampleCount = 160;

        std::ofstream csv(outCsv);
        csv << std::fixed << std::setprecision(10);

        csv
            << "r,"
            << "LambdaL2,"
            << "Eabs,"
            << "VirtualGridResolution,"
            << "DeltaXNorm"
            << "\n";

        for (uint32_t i = 0; i < sampleCount; ++i)
        {
            float t = float(i) / float(sampleCount - 1);

            // Log-spaced r because the metric changes rapidly near the surface.
            float r = rMin * std::pow(rMax / rMin, t);

            float3 q = r * qDir;
            float rInv = 1.0f / r;
            float rInvSq = rInv * rInv;
            float cosXi = -(dot(n, q)) * rInv;

            float3 gradOmega = gradientOmega(q, n, rInv, cosXi, Omega_i);
            float3x3 H_Omega = hessianOmega(q, n, rInv, cosXi, Omega_i);

            std::array<float, 9> ylm;
            std::array<float3, 9> gradYlm;
            std::array<float3x3, 9> hessYlm;

            SHGradientAndHessianL2(qDir, ylm, gradYlm, hessYlm);

            std::array<float3x3, 9> HTotal;

            for (int basisIdx = 0; basisIdx < 9; ++basisIdx)
            {
                HTotal[basisIdx] = float3x3::zeros();

                for (int j = 0; j < 3; ++j)
                {
                    for (int k = j; k < 3; ++k)
                    {
                        float term1 = H_Omega[j][k] * ylm[basisIdx];

                        float term2 = -rInv * (
                            gradOmega[j] * gradYlm[basisIdx][k] +
                            gradOmega[k] * gradYlm[basisIdx][j]
                            );

                        float term3 = Omega_i * rInvSq * hessYlm[basisIdx][j][k];

                        float hessContrib = L * (term1 + term2 + term3);

                        HTotal[basisIdx][j][k] += hessContrib;

                        if (j != k)
                            HTotal[basisIdx][k][j] += hessContrib;
                    }
                }
            }

            float lambdaL2 = lambdaL2OverBands(HTotal);
            float eabs = 0.5f * lambdaL2 * deltaXNorm * deltaXNorm;

            csv
                << r << ","
                << lambdaL2 << ","
                << eabs << ","
                << virtualGridResolution << ","
                << deltaXNorm
                << "\n";
        }

        csv.close();

        logInfo("Wrote analytic error-vs-distance test: " + outCsv);
    }
}

extern "C" FALCOR_API_EXPORT void registerPlugin(Falcor::PluginRegistry& registry)
{
    registry.registerClass<RenderPass, PrecomputeSHCoefficients>();
}

PrecomputeSHCoefficients::PrecomputeSHCoefficients(ref<Device> pDevice, const Properties& props) : RenderPass(pDevice)
{

    for (const auto& [key, value] : props)
    {
        if (key == kShowReconstructedEnvMap)
            mbShowReconstructedEnvMap = value;
        if (key == kShowSHGrid)
            mbShowAdaptiveGrid = value;
        if (key == kShowEdgeAddedVoxels)
            mbShowEdgeAddedVoxels = value;
    }

    mpFbo = Fbo::create(mpDevice);
    Sampler::Desc samplerDesc;
    samplerDesc.setFilterMode(TextureFilteringMode::Linear, TextureFilteringMode::Linear, TextureFilteringMode::Linear);

    mpLinearSampler = mpDevice->createSampler(samplerDesc);
   
}

Properties PrecomputeSHCoefficients::getProperties() const
{
    Properties props;
    props[kShowReconstructedEnvMap] = mbShowReconstructedEnvMap;
    props[kShowSHGrid] = mbShowAdaptiveGrid;
    props[kShowEdgeAddedVoxels] = mbShowEdgeAddedVoxels;
    return props;
}

RenderPassReflection PrecomputeSHCoefficients::reflect(const CompileData& compileData)
{
    // Define the required resources here
    RenderPassReflection reflector;
    const uint2 sz = RenderPassHelpers::calculateIOSize(mOutputSizeSelection, mFixedOutputSize, compileData.defaultTexDims);
    // REMARK MSAA is set via texture sample count. Note that all fbo attachment have to have same sample count.
    reflector.addOutput("output", "Color").texture2D(sz.x, sz.y, 4).format(ResourceFormat::RGBA32Float);
    // Add the required depth output. This always exists.
     reflector.addOutput("depth", "Depth buffer")
         .format(ResourceFormat::D32Float)
         .bindFlags(ResourceBindFlags::DepthStencil)
         .texture2D(sz.x, sz.y, 4);

    //reflector.addOutput("output", "Color").texture2D(sz.x, sz.y).format(ResourceFormat::RGBA16Float);
    //reflector.addOutput("depth", "Depth buffer")
    //    .format(ResourceFormat::D32Float)
    //    .bindFlags(ResourceBindFlags::DepthStencil)
    //    .texture2D(sz.x, sz.y);
    return reflector;
}

void PrecomputeSHCoefficients::execute(RenderContext* pRenderContext, const RenderData& renderData)
{
    auto pTargetFbo = renderData.getTexture("output");
    const float4 clearColor(0, 0, 0, 1);
    mpFbo->attachColorTarget(pTargetFbo, 0);

    // Update frame dimension based on render pass output.
    auto pDepth = renderData.getTexture("depth");

    //  Clear depth buffer.
    //pRenderContext->clearDsv(pDepth->getDSV().get(), 1.f, 0);
    mpFbo->attachDepthStencilTarget(pDepth);

    pRenderContext->clearFbo(mpFbo.get(), clearColor, 1.0f, 0, FboAttachmentType::Color);

    if (mpScene)
    {
#pragma region  ====Finite difference verification code block (Row-by-Row)===
//if (mbVerify)
//{
//    // 1. Open CSV File
//    std::ofstream csvFile("verification_data.csv");
//    // Added AnalyticHessXZ and FDHessXZ columns
//    csvFile << "WorldX,WorldZ,AnalyticGrad,FDGrad,AnalyticHessXX,FDHessXX,AnalyticHessXZ,FDHessXZ\n";
//
//    // 2. Settings
//    int basisIdx = 4; // l=2, m=-2 (xy in world, depends on basis def)
//    int numBasis = 9;
//    float verificationHSq = verificationH * verificationH;
//    float verificationHSq4 = 4.0f * verificationHSq; // Denominator for Mixed FD (4h^2)
//
//    // Step size for the grid coverage
//    double step = (2.0 * verificationExtent) / (double)(verificationRes - 1);
//
//    // =================================================================================
//    // MAIN LOOP: Process Row by Row (Z-Axis)
//    // =================================================================================
//    for (uint32_t zIdx = 0; zIdx < verificationRes; ++zIdx)
//    {
//        // A. Generate ONLY CENTER Positions (World Space)
//        std::vector<float3> rowPositions;
//        rowPositions.clear();
//        rowPositions.reserve(verificationRes);
//
//        double z = -verificationExtent + step * (double)zIdx;
//
//        for (uint32_t xIdx = 0; xIdx < verificationRes; ++xIdx)
//        {
//            double x = -verificationExtent + step * (double)xIdx;
//            rowPositions.push_back(float3((float)x, verificationY, (float)z));
//        }
//
//        int numProbesInRow = (int)rowPositions.size();
//
//        // B. Update GPU Buffers
//        mpProbePosBuffer = mpDevice->createStructuredBuffer(
//            sizeof(float3), numProbesInRow,
//            ResourceBindFlags::ShaderResource, MemoryType::DeviceLocal,
//            rowPositions.data()
//        );
//
//        mpProbeSamplingResultBuffer = mpDevice->createStructuredBuffer(
//            sizeof(ProbeSampleData),
//            numSamplesPerProbe * numProbesInRow,
//            ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess,
//            MemoryType::DeviceLocal
//        );
//
//        // C. Dispatch Ray Tracing
//        auto rtVar = mpRtVars->getRootVar();
//        rtVar["gProbeDirSamples"] = mpProbeDirSamplesBuffer;
//        rtVar["gProbePositions"] = mpProbePosBuffer;
//        if (mpEmissiveSampler) mpEmissiveSampler->bindShaderData(rtVar["PerFrameCB"]["emissiveSampler"]);
//
//        rtVar["gProbeSamplingOutput"] = mpProbeSamplingResultBuffer;
//        rtVar["PerFrameCB"]["numSamplePerProbe"] = numSamplesPerProbe;
//        rtVar["PerFrameCB"]["sampleIndex"] = mSampleIndex++;
//
//        mpScene->raytrace(pRenderContext, mpRtProgram.get(), mpRtVars, uint3(numSamplesPerProbe, numProbesInRow, 1));
//
//        // D. Readback Row Data
//        ProbeSampleData* rowAllData = new ProbeSampleData[numSamplesPerProbe * numProbesInRow];
//        mpProbeSamplingResultBuffer->getBlob(rowAllData, 0, numSamplesPerProbe * numProbesInRow * sizeof(ProbeSampleData));
//
//        // E. Process Probes
//        for (int i = 0; i < verificationRes; ++i)
//        {
//            // 1. Extract Samples for the CENTER Probe
//            std::vector<ProbeSampleData> centerSamples;
//            centerSamples.reserve(numSamplesPerProbe);
//            int offset = i * numSamplesPerProbe;
//            for (int s = 0; s < numSamplesPerProbe; ++s) centerSamples.push_back(rowAllData[offset + s]);
//
//            // 2. Define Center Position (World Space)
//            float3 centerPosWorld = rowPositions[i];
//
//            // 3. Define Polar Coordinates (Swizzled: Z, X, Y)
//            float3 centerPosPolar;
//            centerPosPolar.x = centerPosWorld.z;
//            centerPosPolar.y = centerPosWorld.x;
//            centerPosPolar.z = centerPosWorld.y;
//
//            // -----------------------------------------------------------------
//            // 4. Stencil Generation (Using Reprojection / Virtual shifts)
//            // -----------------------------------------------------------------
//
//            // --- A. Pure X Derivative Stencil (For Gradient and HessXX) ---
//            // World X maps to Polar Y
//            float3 posPolarXPlus = centerPosPolar; posPolarXPlus.y += verificationH;
//            float3 posPolarXMinus = centerPosPolar; posPolarXMinus.y -= verificationH;
//
//            // --- B. Mixed XZ Derivative Stencil (For HessXZ) ---
//            // World X maps to Polar Y, World Z maps to Polar X
//            float3 posPolarXPZP = centerPosPolar; // X+ Z+
//            posPolarXPZP.y += verificationH; posPolarXPZP.x += verificationH;
//
//            float3 posPolarXPZM = centerPosPolar; // X+ Z-
//            posPolarXPZM.y += verificationH; posPolarXPZM.x -= verificationH;
//
//            float3 posPolarXMZP = centerPosPolar; // X- Z+
//            posPolarXMZP.y -= verificationH; posPolarXMZP.x += verificationH;
//
//            float3 posPolarXMZM = centerPosPolar; // X- Z-
//            posPolarXMZM.y -= verificationH; posPolarXMZM.x -= verificationH;
//
//
//            // 5. Compute Coefficients (Numerical Reprojection)
//            // Note: We use 'centerSamples' for ALL calls to lock the noise.
//
//            float coeffR = calculateChannelRSHCoeffLM(basisIdx, centerSamples);
//
//            // Pure X
//            float coeffXP = calculateCoeffRPrime(centerSamples, centerPosPolar, posPolarXPlus, numBasis, basisIdx);
//            float coeffXM = calculateCoeffRPrime(centerSamples, centerPosPolar, posPolarXMinus, numBasis, basisIdx);
//
//            // Mixed XZ
//            float coeffXPZP = calculateCoeffRPrime(centerSamples, centerPosPolar, posPolarXPZP, numBasis, basisIdx);
//            float coeffXPZM = calculateCoeffRPrime(centerSamples, centerPosPolar, posPolarXPZM, numBasis, basisIdx);
//            float coeffXMZP = calculateCoeffRPrime(centerSamples, centerPosPolar, posPolarXMZP, numBasis, basisIdx);
//            float coeffXMZM = calculateCoeffRPrime(centerSamples, centerPosPolar, posPolarXMZM, numBasis, basisIdx);
//
//
//            // 6. Finite Difference Calculations
//
//            // Gradient d/dX (World)
//            float fdGrad = (coeffXP - coeffR) / verificationH; // Forward diff 
//
//            // Hessian d^2/dX^2 (World)
//            float fdHessXX = (coeffXP - 2.0f * coeffR + coeffXM) / verificationHSq;
//
//            // Hessian d^2/dXdZ (World) -> Formula: (f_++ - f_+- - f_-+ + f_--) / 4h^2
//            float fdHessXZ = (coeffXPZP - coeffXPZM - coeffXMZP + coeffXMZM) / verificationHSq4;
//
//
//            // 7. Analytic Calculations (Ours)
//            float3 analyticGrad = float3(0, 0, 0);
//            float3x3 analyticHess = float3x3::zeros();
//            calculateChannelRGradAndHessianSHCoeffLM(centerPosPolar, centerSamples, samplingDirs, basisIdx, analyticGrad, analyticHess);
//
//            // Mapping:
//            // AnalyticGrad.y corresponds to d/dPolarY -> d/dWorldX
//            // AnalyticHess[1][1] corresponds to d^2/dPolarY^2 -> d^2/dWorldX^2
//            // AnalyticHess[0][1] corresponds to d^2/dPolarX dPolarY -> d^2/dWorldZ dWorldX
//
//            // 8. Write to CSV
//            csvFile << centerPosWorld.x << "," << centerPosWorld.z << ","
//                << analyticGrad.y << "," << fdGrad << ","
//                << analyticHess[1][1] << "," << fdHessXX << ","
//                << analyticHess[0][1] << "," << fdHessXZ << "\n";
//        }
//
//        delete[] rowAllData; // Cleanup array
//    }
//
//    csvFile.close();
//    mbVerify = false;
//    mSampleIndex++;
//    std::cout << "Verification Data Export Complete!" << std::endl;
//}
#pragma endregion
#pragma region  ====Uniform grid SH coeff block===
        //if (!mbFinishSHPrecompute)
        //{
        //    auto rtVar = mpRtVars->getRootVar();
        //    rtVar["gProbeDirSamples"] = mpProbeDirSamplesBuffer;
        //    rtVar["gProbePositions"] = mpProbePosBuffer;

        //    if (mpEmissiveSampler)
        //        mpEmissiveSampler->bindShaderData(rtVar["PerFrameCB"]["emissiveSampler"]);

        //    rtVar["gProbeSamplingOutput"] = mpProbeSamplingResultBuffer;
        //    rtVar["PerFrameCB"]["numSamplePerProbe"] = numSamplePerProbe;
        //    rtVar["PerFrameCB"]["sampleIndex"] = mSampleIndex++;
        //    int numProbe = mProbeGrid.resolution.x * mProbeGrid.resolution.y * mProbeGrid.resolution.z;
        //    mpScene->raytrace(pRenderContext, mpRtProgram.get(), mpRtVars, uint3(numSamplePerProbe, numProbe, 1));

        //    ProbeSampleData* samplingData = new ProbeSampleData[numSamplePerProbe * numProbe];
        //    mpProbeSamplingResultBuffer->getBlob(samplingData, 0, numSamplePerProbe * numProbe * sizeof(ProbeSampleData));

        //    //// Now pData points to your results, size is sampleCount
        //    //for (int i = 0; i < sampleCount; ++i)
        //    //{
        //    //    float4 result = pData[i];
        //    //    // Process result as needed
        //    //    logInfo(fmt::format("AAAA: {:.3f} {:.3f} {:.3f}", result.x, result.y, result.z));
        //    //}

        //    mProbeGrid.probesSHCoeffs.clear();
        //    mProbeGrid.probesSHCoeffs.reserve(mProbeGrid.numBasis*numProbe);

        //    mProbeGrid.probesSHCoeffsGradients.clear();
        //    mProbeGrid.probesSHCoeffsGradients.reserve(mProbeGrid.numBasis * numProbe);


        //    mProbeGrid.probesSHCoeffsHessians.clear();
        //    mProbeGrid.probesSHCoeffsHessians.reserve(mProbeGrid.numBasis * numProbe);

        //    for (int probeIdx = 0; probeIdx < numProbe; ++probeIdx)
        //    {
        //        int offset = probeIdx * numSamplePerProbe;
        //        std::vector<float4> shCoeffs;
        //        std::vector<ProbeSampleData> probeSamplingResults;
        //        probeSamplingResults.clear();
        //        probeSamplingResults.reserve(numSamplePerProbe);
        //        for (int i = 0; i < numSamplePerProbe; i++)
        //        {
        //            probeSamplingResults.push_back(samplingData[offset + i]);
        //        }
        //        calculateSHCoeffs(shCoeffs, probeSamplingResults, numSamplePerProbe);

        //        std::vector<GradSHCoeff> shCoeffsGradients;
        //        std::vector<HessianSHCoeff> shCoeffsHessians;
        //        float3 xPolar;
        //        xPolar.x = mProbeGrid.probesPos[probeIdx].z;
        //        xPolar.y = mProbeGrid.probesPos[probeIdx].x;
        //        xPolar.z = mProbeGrid.probesPos[probeIdx].y;
        //        calculateSHCoeffsGradientsAndHessians(shCoeffsGradients, shCoeffsHessians, xPolar, probeSamplingResults);
        //        mProbeGrid.probesSHCoeffs.insert(mProbeGrid.probesSHCoeffs.end(), shCoeffs.begin(), shCoeffs.end());
        //        mProbeGrid.probesSHCoeffsGradients.insert(mProbeGrid.probesSHCoeffsGradients.end(), shCoeffsGradients.begin(), shCoeffsGradients.end());
        //        mProbeGrid.probesSHCoeffsHessians.insert(mProbeGrid.probesSHCoeffsHessians.end(), shCoeffsHessians.begin(), shCoeffsHessians.end());
        //    }

        //    //if (sceneName == "arcade")
        //    //{
        //    //    saveProbeGridToFile(mProbeGrid, "ProbeGridArcade.txt");
        //    //}
        //    //else
        //    //{
        //    //     saveProbeGridToFile(mProbeGrid, "ProbeGridCornell.txt");
        //    //}

        //    if (sceneName == "arcade")
        //    {
        //        saveProbeGridToFileWithGradAndHessian(mProbeGrid, "ProbeGridWithGradAndHessianArcade.txt");
        //    }
        //    else
        //    {
        //        saveProbeGridToFileWithGradAndHessian(mProbeGrid, "ProbeGridWithGradAndHessianCornell.txt");
        //    }

        //    mpGridSHCoeffsBuffer = mpDevice->createStructuredBuffer(
        //        sizeof(float4),
        //        mProbeGrid.numBasis * numProbe,
        //        ResourceBindFlags::ShaderResource,
        //        MemoryType::DeviceLocal,
        //        mProbeGrid.probesSHCoeffs.data()
        //    );
        //    mpGridSHCoeffsBuffer->setName("SH Grid Coeffs");

        //    mbFinishSHPrecompute = true;
        //    delete[] samplingData;
        //}
#pragma endregion
#pragma region  ====Adaptive probe volume construction block===
        // build adaptive probe grid here
#if CURRENT_PROBE_MODE == PROBE_MODE_ADAPTIVE
        if (mNeedRebuildProbeVolume)
        {
            using clock = std::chrono::high_resolution_clock;
            auto tStart = clock::now();
           
            //ProgressiveRefineBuild(pRenderContext);
            SinglePassBuild(pRenderContext);
            auto tEnd = clock::now();
            double ms = std::chrono::duration<double, std::milli>(tEnd - tStart).count();
            mAdaptiveProbeVolume->setBuildTimeMs(ms);
            // 3. TODO: upload to GPU for visualization
            mAdaptiveProbeVolume->uploadToGPU();

            //std::string debugFileName = saveToFileName;
            //size_t pos = debugFileName.find_last_of('.');
            //debugFileName.replace(pos, std::string::npos, "_DebugViz.txt");
            std::string debugFileName =
                replaceExtensionOrAppend(
                    saveToFileName,
                    "_DebugViz.txt"
                );

            std::string residualLogFileName =
                replaceExtensionOrAppend(
                    saveToFileName,
                    "_ResidualStats.txt"
                );

            std::string residualSummaryCsvFileName =
                replaceExtensionOrAppend(
                    saveToFileName,
                    "_ResidualSummary.csv"
                );


            std::string residualLevelCsvFileName =
                replaceExtensionOrAppend(
                    saveToFileName,
                    "_ResidualByLevel.csv"
                );

            // Separate residual statistics log.
            mAdaptiveProbeVolume->writeResidualPaperStatsLog(
                residualLogFileName
            );

            std::string residualMainStatsFileName =
                replaceExtensionOrAppend(
                    saveToFileName,
                    "_ResidualMainStats.txt"
                );

            std::string residualMainStatsCsvFileName =
                replaceExtensionOrAppend(
                    saveToFileName,
                    "_ResidualMainStats.csv"
                );

            mAdaptiveProbeVolume->writeResidualMainStatsLog(
                residualMainStatsFileName
            );

            mAdaptiveProbeVolume->exportResidualMainStatsCSV(
                residualMainStatsCsvFileName
            );
            // Separate residual CSV outputs.
            mAdaptiveProbeVolume->exportResidualPaperStatsCSV(
                residualSummaryCsvFileName,
                residualLevelCsvFileName
            );


            logInfo("Wrote residual statistics log: " + residualLogFileName);
            logInfo("Wrote residual summary CSV: " + residualSummaryCsvFileName);
            logInfo("Wrote residual by-level CSV: " + residualLevelCsvFileName);

            mAdaptiveProbeVolume->printDebugInfo(debugFileName);
            logInfo("Wrote hierarchy debug file: " + debugFileName);
            mAdaptiveProbeVolume->saveToFile(saveToFileName);
            mNeedRebuildProbeVolume = false;
            mpProbeVisualizePass->setVolumeData(mAdaptiveProbeVolume->getProbes());

        }
#endif
#pragma endregion
#pragma region  ====Uniform probe volume construction block===
#if CURRENT_PROBE_MODE == PROBE_MODE_UNIFORM
        if (mNeedRebuildProbeVolume)
        {
            using clock = std::chrono::high_resolution_clock;
            auto tStart = clock::now();
            // 1. Initialize Grid Structure (Calculates resolution and total probes)
            mUniformProbeVolume->initGrid(mpScene, unifromGridSize);
            uint3 probeCountDim = mUniformProbeVolume->getProbeCountDim(); // (N+1) corners

            // 2. Prepare reusable buffers for a single Row (X-dimension)
            uint32_t numProbesPerRow = probeCountDim.x;

            // Position Buffer for one row
            mpProbePosBuffer = mpDevice->createStructuredBuffer(
                sizeof(float3), numProbesPerRow, ResourceBindFlags::ShaderResource, MemoryType::DeviceLocal
            );

            // Sampling result buffer for one row
            mpProbeSamplingResultBuffer = mpDevice->createStructuredBuffer(
                sizeof(ProbeSampleData),
                numSamplesPerProbe * numProbesPerRow,
                ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess,
                MemoryType::DeviceLocal
            );

            // Temp storage for CPU readback of one row
            std::vector<ProbeSampleData> rowSamplingData(numSamplesPerProbe * numProbesPerRow);
            std::vector<float3> rowPositions(numProbesPerRow);

            // --------------------------------------------------------------------------
            // 3. Row-by-Row Processing Loop
            // --------------------------------------------------------------------------
            for (uint32_t z = 0; z < probeCountDim.z; ++z)
            {
                for (uint32_t y = 0; y < probeCountDim.y; ++y)
                {
                    // A. Calculate positions for the current row (X-axis)
                    for (uint32_t x = 0; x < probeCountDim.x; ++x)
                    {
                        // Position = Origin + Index * CellSize
                        rowPositions[x] = mUniformProbeVolume->getMinPoint() +
                            float3(x, y, z) * mUniformProbeVolume->getCellSize();
                    }

                    // B. Upload row positions to GPU
                    mpProbePosBuffer->setBlob(rowPositions.data(), 0, numProbesPerRow * sizeof(float3));

                    // C. Dispatch Ray Tracing for this row
                    auto rtVar = mpRtVars->getRootVar();
                    rtVar["gProbeDirSamples"] = mpProbeDirSamplesBuffer;
                    rtVar["gProbePositions"] = mpProbePosBuffer;
                    rtVar["gProbeSamplingOutput"] = mpProbeSamplingResultBuffer;
                    rtVar["PerFrameCB"]["sampleIndex"] = mSampleIndex++;
                    rtVar["PerFrameCB"]["numSamplePerProbe"] = numSamplesPerProbe;

                    if (mpEmissiveSampler)
                        mpEmissiveSampler->bindShaderData(rtVar["PerFrameCB"]["emissiveSampler"]);

                    mpScene->raytrace(pRenderContext, mpRtProgram.get(), mpRtVars, uint3(numSamplesPerProbe, numProbesPerRow, 1));
                    pRenderContext->submit(true);
                    // D. Synchronize and Readback row results
                    mpProbeSamplingResultBuffer->getBlob(rowSamplingData.data(), 0, rowSamplingData.size() * sizeof(ProbeSampleData));

                    // E. Process each probe in the row on CPU
                    for (uint32_t x = 0; x < numProbesPerRow; ++x)
                    {
                        uint32_t probeIdx = (z * probeCountDim.y * probeCountDim.x) + (y * probeCountDim.x) + x;

                        //int offset = x * numSamplesPerProbe;
                        //std::vector<ProbeSampleData> probeSamplingResults;
                        //probeSamplingResults.assign(rowSamplingData.begin() + offset, rowSamplingData.begin() + offset + numSamplesPerProbe);

                        // Convert to Polar for SH Math (Z, X, Y mapping)
                        float3 xPolar = float3(rowPositions[x].z, rowPositions[x].x, rowPositions[x].y);

                        std::vector<float3> coeffs;
                        std::vector<GradSHCoeff> grads;
                        size_t offset = size_t(x) * size_t(numSamplesPerProbe);
                        const ProbeSampleData* samples =
                            rowSamplingData.data() + offset;

                        // Perform Physics Calculations
                        calculateSHCoeffsGradients(grads, xPolar, samples, numSamplesPerProbe, samplingDirs);
                        calculateSHCoeffs(coeffs, samples, numSamplesPerProbe);

                        // Store in the Volume object
                        mUniformProbeVolume->setProbeData(probeIdx, coeffs, grads);
                    }
                }
                mpDevice->wait();
            }
            auto tEnd = clock::now();
            double ms = std::chrono::duration<double, std::milli>(tEnd - tStart).count();
            mUniformProbeVolume->setBuildTimeMs(ms);
            // --------------------------------------------------------------------------
            // 4. Finalize
            // --------------------------------------------------------------------------
            mUniformProbeVolume->uploadToGPU();
            mUniformProbeVolume->saveToFile(saveToFileName);

            mpProbeVisualizePass->setUniformVolumeData(
                mUniformProbeVolume->getMinPoint(),
                mUniformProbeVolume->getCellSize(),
                mUniformProbeVolume->getCellResolution()
            );

            mNeedRebuildProbeVolume = false;
        }
#endif
#pragma endregion
        // visualize probes
        //if (mbFinishSHPrecompute)
        //{
            pRenderContext->clearDsv(pDepth->getDSV().get(), 1.f, 0);

            auto shShaderRootVar = mpVars->getRootVar();
             shShaderRootVar["gLinearSampler"] = mpLinearSampler;
 #if CURRENT_PROBE_MODE == PROBE_MODE_ADAPTIVE
             shShaderRootVar["gCornerBuffer"] = mAdaptiveProbeVolume->getCornerBuffer();
             shShaderRootVar["gProbeBuffer"] = mAdaptiveProbeVolume->getProbeBuffer();
#else
             mUniformProbeVolume->bindShaderData(shShaderRootVar);
#endif
             mpScene->rasterize(pRenderContext, mpGraphicsState.get(), mpVars.get(), mpRasterState, mpRasterState);

             float4x4 viewProj = mpScene->getCamera()->getViewProjMatrix();

             if (mbShowAdaptiveGrid && mpProbeVisualizePass)
             {
                 mpProbeVisualizePass->setCameraData(viewProj);
                 mpProbeVisualizePass->execute(pRenderContext, mpFbo);
             }
        //}
    }
}
float PrecomputeSHCoefficients::calculateCoeffRPrime(std::vector<ProbeSampleData> probeSamplingResults, float3 xPolar, float3 xPolarXPrime, int numBasis, int basisIdx)
{
    // refer to Krivanek 2005 sect 4.3.1
    std::vector<float> XPrimeSHBasisTable(9 * numSamplesPerProbe);
    std::vector<float> XPrimeSolidAngleWeight(numSamplesPerProbe);
    std::vector<float3> XPrimeSHGradientTable(9 * numSamplesPerProbe);
    std::vector<float3x3> XPrimeSHHessianTable(9 * numSamplesPerProbe);

    XPrimeSHBasisTable.assign(9 * numSamplesPerProbe, 0.0f);
    XPrimeSolidAngleWeight.assign(numSamplesPerProbe, 0.0f);
    XPrimeSHGradientTable.assign(9 * numSamplesPerProbe, float3(0.0f));
    XPrimeSHHessianTable.assign(9 * numSamplesPerProbe, float3x3::zeros());

    std::vector<float3> xPrimeDirs;
    xPrimeDirs.clear();
    xPrimeDirs.reserve(numSamplesPerProbe);

    // Solid angle (uniform per patch)
    float Omega_i = (4.0f * (float)M_PI) / (float)numSamplesPerProbe;

    // calculate new directions from x' to all sample points
    for (int sampleIdx = 0; sampleIdx < numSamplesPerProbe; ++sampleIdx)
    {
        float3 sPolar =
            float3(probeSamplingResults[sampleIdx].s.x, probeSamplingResults[sampleIdx].s.y, probeSamplingResults[sampleIdx].s.z);
        float3 nPolar =
            float3(probeSamplingResults[sampleIdx].n.x, probeSamplingResults[sampleIdx].n.y, probeSamplingResults[sampleIdx].n.z);
        float3 L = float3(probeSamplingResults[sampleIdx].Li.x, probeSamplingResults[sampleIdx].Li.y, probeSamplingResults[sampleIdx].Li.z);

        float3 xPrimeNormDir = normalize(sPolar - xPolarXPrime); // direction from x' to same sample point
        xPrimeDirs.push_back(xPrimeNormDir);

        if (probeSamplingResults[sampleIdx].hitT < 0.0f) // ray miss
        {
            XPrimeSolidAngleWeight[sampleIdx] = 0.0f;
        }
        else
        {
            // Geometry and direction
            float3 qk = sPolar - xPolar;
            float rk = math::length(qk);

            // cosξ = -(n · q) / r
            float cosXik = -(math::dot(nPolar, qk)) / rk;

            float3 qkPrime = sPolar - xPolarXPrime;
            float rkPrime = math::length(qkPrime);
            float cosXikPrime = -(math::dot(nPolar, qkPrime)) / rkPrime;

            XPrimeSolidAngleWeight[sampleIdx] = (rk * rk * cosXikPrime) / (rkPrime * rkPrime * cosXik);
        }
    }

     // calculate SH basis and derivatives at x' using new directions
    initSHBasisGradientAndHessianTables(xPrimeDirs, XPrimeSHBasisTable, XPrimeSHGradientTable, XPrimeSHHessianTable);

    // calculate SH coeff at x' using new directions and solid angle correction
    float coeffRPrime = 0.0f;
    // For each direction sample
    for (int sampleIdx = 0; sampleIdx < numSamplesPerProbe; ++sampleIdx)
    {
        const float4& sample = probeSamplingResults[sampleIdx].Li;
        float shYPrime = XPrimeSHBasisTable[numBasis * sampleIdx + basisIdx]; // SH basis value for this direction

        coeffRPrime += sample.r * shYPrime * XPrimeSolidAngleWeight[sampleIdx];
    }
    coeffRPrime *= Omega_i;

    return coeffRPrime;
}

void PrecomputeSHCoefficients::SinglePassBuild(RenderContext* pRenderContext)
{
    // 1. Initialize
           //mAdaptiveProbeVolume->startBuild(mpScene, ErrorThreshold, useRelativeError);

    const uint3 seedResolution = uint3(1, 1, 1);
    //const uint3 seedResolution = uint3(4,4,4); 
    //const uint3 seedResolution = uint3(8,8,8); 
    //const uint3 seedResolution = uint3(16, 16, 16);
    //mAdaptiveProbeVolume->setErrorMetricMode(
    //    useIrradianceSpaceMetric
    //    ? AdaptiveProbeVolume::ErrorMetricMode::IrradianceSpace
    //    : AdaptiveProbeVolume::ErrorMetricMode::SHSpace
    //);


    mAdaptiveProbeVolume->setResidualCorrection(
        useResidualCorrection,
        residualPruneStrength,
        residualRefineStrength,
        residualCorrectionMinScale,
        residualCorrectionMaxScale,
        residualConfidenceTau0,
        residualConfidenceTau1,
        residualConfidenceEps
    );

    mAdaptiveProbeVolume->startBuildSeeded(mpScene, seedResolution, ErrorThreshold, useRelativeError);

    //const uint32_t kMaxCornersPerDispatch = 8192; // tune this
    const uint32_t kMaxCornersPerDispatch = 4096; // tune this
    //const uint32_t kMaxCornersPerDispatch = 1024; // tune this

    while (mAdaptiveProbeVolume->hasPendingBatch())
    {
        uint32_t totalPending = mAdaptiveProbeVolume->getPendingCornerCount();

        for (uint32_t batchStart = 0; batchStart < totalPending; batchStart += kMaxCornersPerDispatch)
        {
            uint32_t batchCount = std::min(kMaxCornersPerDispatch, totalPending - batchStart);

            std::vector<float3> pendingProbePositions;
            mAdaptiveProbeVolume->getPendingPositionsRange(batchStart, batchCount, pendingProbePositions);

            uint32_t numProbes = (uint32_t)pendingProbePositions.size();

            mpProbePosBuffer = mpDevice->createStructuredBuffer(
                sizeof(float3), numProbes, ResourceBindFlags::ShaderResource, MemoryType::DeviceLocal, pendingProbePositions.data()
            );

            mpProbeSamplingResultBuffer = mpDevice->createStructuredBuffer(
                sizeof(ProbeSampleData),
                numSamplesPerProbe * numProbes,
                ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess,
                MemoryType::DeviceLocal
            );

            auto rtVar = mpRtVars->getRootVar();
            rtVar["gProbeDirSamples"] = mpProbeDirSamplesBuffer;
            rtVar["gProbePositions"] = mpProbePosBuffer;
            rtVar["PerFrameCB"]["sampleIndex"] = mSampleIndex++;
            if (mpEmissiveSampler)
                mpEmissiveSampler->bindShaderData(rtVar["PerFrameCB"]["emissiveSampler"]);

            rtVar["gProbeSamplingOutput"] = mpProbeSamplingResultBuffer;
            rtVar["PerFrameCB"]["numSamplePerProbe"] = numSamplesPerProbe;

            mpScene->raytrace(pRenderContext, mpRtProgram.get(), mpRtVars, uint3(numSamplesPerProbe, numProbes, 1));

            std::vector<ProbeSampleData> allProbeSamplingData(numSamplesPerProbe * numProbes);
            mpProbeSamplingResultBuffer->getBlob(
                allProbeSamplingData.data(),
                0,
                numSamplesPerProbe * numProbes * sizeof(ProbeSampleData)
            );

            std::vector<std::vector<float3>> coeffsBatch(numProbes);
            std::vector<std::vector<GradSHCoeff>> gradsBatch(numProbes);
            std::vector<std::vector<float3x3>> hessiansBatch(numProbes);

            for (uint32_t probeIdx = 0; probeIdx < numProbes; ++probeIdx)
            {
                int offset = probeIdx * numSamplesPerProbe;
                //std::vector<ProbeSampleData> probeSamplingResults;
                //probeSamplingResults.reserve(numSamplesPerProbe);
                //for (int sampleIdx = 0; sampleIdx < numSamplesPerProbe; ++sampleIdx)
                //    probeSamplingResults.push_back(allProbeSamplingData[offset + sampleIdx]);

                const ProbeSampleData* samples =
                    allProbeSamplingData.data() + offset;
                float3 xPolar;
                xPolar.x = pendingProbePositions[probeIdx].z;
                xPolar.y = pendingProbePositions[probeIdx].x;
                xPolar.z = pendingProbePositions[probeIdx].y;

                calculateSHCoeffsGradientsRGBAndHessiansLum(
                    gradsBatch[probeIdx],
                    hessiansBatch[probeIdx],
                    xPolar,
                    samples,
                    numSamplesPerProbe,
                    samplingDirs
                );
                calculateSHCoeffs(coeffsBatch[probeIdx], samples, numSamplesPerProbe);
            }

            mAdaptiveProbeVolume->setCornerDataRange(batchStart, coeffsBatch, gradsBatch, hessiansBatch);
        }
        mAdaptiveProbeVolume->finishBatch();
    }
}

void PrecomputeSHCoefficients::ProgressiveRefineBuild(RenderContext* pRenderContext)
{
    //const uint32_t kMaxCornersPerDispatch = 1024;
    //const uint32_t kMaxCornersPerDispatch = 4096;
    const uint32_t kMaxCornersPerDispatch = 8192;

    struct LocalProgressiveStage
    {
        uint32_t spp;
        float thresholdScale;
        bool storeFullRuntimeData;
        bool refineAfterThisStage;
    };

    const LocalProgressiveStage stages[] =
    {
        {128,  1.0f, false, true}, // coarse metric-only topology
        {1024, 1.0f, true,  false}  // full adaptive: store runtime + refine onward
    };

    mAdaptiveProbeVolume->resetBuildStats();

    const uint3 seedResolution = uint3(1, 1, 1);
    mAdaptiveProbeVolume->startBuildSeeded(mpScene, seedResolution, ErrorThreshold, useRelativeError);

    // Reusable buffers.
    mpProbePosBuffer = mpDevice->createStructuredBuffer(
        sizeof(float3),
        kMaxCornersPerDispatch,
        ResourceBindFlags::ShaderResource,
        MemoryType::DeviceLocal
    );

    mpProbeSamplingResultBuffer = mpDevice->createStructuredBuffer(
        sizeof(ProbeSampleData),
        kMaxSamplesPerProbe * kMaxCornersPerDispatch,
        ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess,
        MemoryType::DeviceLocal
    );

    std::vector<ProbeSampleData> allProbeSamplingData;
    allProbeSamplingData.resize(uint64_t(kMaxSamplesPerProbe) * kMaxCornersPerDispatch);

    for (int stageIdx = 0; stageIdx < 2; ++stageIdx)
    {
        const auto& stage = stages[stageIdx];

        // ------------------------------------------------------------
        // Leaf-only rule:
        // Stage 0 starts from root pending corners.
        // Later stages recheck current leaves only.
        // ------------------------------------------------------------
        //if (stageIdx > 0)
        //{
        //    if (stage.storeFullRuntimeData)
        //        mAdaptiveProbeVolume->scheduleLeafCornersForRefinementRecheck();
        //    else
        //        mAdaptiveProbeVolume->scheduleLeafCornersForRefinementRecheck();
        //}

        if (stageIdx == 1)
        {
            mAdaptiveProbeVolume->scheduleAllLeafCornersForRuntimeBake();
        }

        mAdaptiveProbeVolume->setCurrentThreshold(ErrorThreshold * stage.thresholdScale);

        while (mAdaptiveProbeVolume->hasPendingBatch())
        {
            uint32_t totalPending = mAdaptiveProbeVolume->getPendingCornerCount();

            for (uint32_t batchStart = 0; batchStart < totalPending; batchStart += kMaxCornersPerDispatch)
            {
                uint32_t batchCount = std::min(kMaxCornersPerDispatch, totalPending - batchStart);

                std::vector<float3> pendingPositions;
                mAdaptiveProbeVolume->getPendingPositionsRange(batchStart, batchCount, pendingPositions);

                uint32_t numCorners = (uint32_t)pendingPositions.size();
                if (numCorners == 0) continue;

                mpProbePosBuffer->setBlob(
                    pendingPositions.data(),
                    0,
                    numCorners * sizeof(float3)
                );

                auto rtVar = mpRtVars->getRootVar();
                rtVar["gProbeDirSamples"] = mpProbeDirSamplesBuffer;
                rtVar["gProbePositions"] = mpProbePosBuffer;
                rtVar["gProbeSamplingOutput"] = mpProbeSamplingResultBuffer;

                rtVar["PerFrameCB"]["numSamplePerProbe"] = stage.spp;
                rtVar["PerFrameCB"]["sampleIndex"] = mSampleIndex++;

                if (mpEmissiveSampler)
                    mpEmissiveSampler->bindShaderData(rtVar["PerFrameCB"]["emissiveSampler"]);

                mpScene->raytrace(
                    pRenderContext,
                    mpRtProgram.get(),
                    mpRtVars,
                    uint3(stage.spp, numCorners, 1)
                );

                mAdaptiveProbeVolume->recordTraceBatch(numCorners, stage.spp);

                mpProbeSamplingResultBuffer->getBlob(
                    allProbeSamplingData.data(),
                    0,
                    uint64_t(stage.spp) * numCorners * sizeof(ProbeSampleData)
                );

                // ------------------------------------------------------------
                // CPU postprocess.
                // Keep vector version for now so it matches your existing SH funcs.
                // ------------------------------------------------------------
                for (uint32_t cornerLocalIdx = 0; cornerLocalIdx < numCorners; ++cornerLocalIdx)
                {
                    uint32_t offset = cornerLocalIdx * stage.spp;

                    const ProbeSampleData* samples =
                        allProbeSamplingData.data() + offset;

                    float3 xPolar;
                    xPolar.x = pendingPositions[cornerLocalIdx].z;
                    xPolar.y = pendingPositions[cornerLocalIdx].x;
                    xPolar.z = pendingPositions[cornerLocalIdx].y;

                    if (stage.storeFullRuntimeData)
                    {
                        std::vector<float3> coeffs;
                        std::vector<GradSHCoeff> grads;
                        std::vector<float3x3> hessians;
                        calculateSHCoeffs(coeffs, samples, stage.spp);

                        //at this point grid is decided only retrace for matching spp with uniform grid
                        calculateSHCoeffsGradients(grads, xPolar, samples, stage.spp, samplingDirs);
                        mAdaptiveProbeVolume->setCornerRuntimeData(
                            batchStart + cornerLocalIdx,
                            coeffs,
                            grads
                        );

                        //calculateSHCoeffsGradientsRGBAndHessiansLum(
                        //    grads,
                        //    hessians,
                        //    xPolar,
                        //    samples,
                        //    stage.spp,
                        //    samplingDirs
                        //);

                         //   mAdaptiveProbeVolume->setCornerData(
                         //   batchStart + cornerLocalIdx,
                         //   coeffs,
                         //   grads,
                         //   hessians
                         //);
                    }
                    else
                    {
                        float coeffVecL2 = 0.0f;
                        float maxLambdaVecL2 = 0.0f;
                        bool isValid = true;

                        calculateSHBuildMetricsOnly(
                            coeffVecL2,
                            maxLambdaVecL2,
                            xPolar,
                            samples,
                            stage.spp,
                            samplingDirs,
                            useRelativeError
                        );

                        mAdaptiveProbeVolume->setCornerMetricData(
                            batchStart + cornerLocalIdx,
                            coeffVecL2,
                            maxLambdaVecL2
                        );
                    }
                }
            }
            if (stage.refineAfterThisStage)
            {
                mAdaptiveProbeVolume->finishBatch();
            }
            else
            {
                mAdaptiveProbeVolume->clearPendingBatch();
            }
        }
        //if (stageIdx == 0)
        //{
        //    mAdaptiveProbeVolume->printCoarseStageDebugInfo("CoarseStageDebug.txt");
        //}
    }
}

void PrecomputeSHCoefficients::renderUI(
    Gui::Widgets& widget
)
{
    if (widget.checkbox(
        "Show Reconstructed Env Map",
        mbShowReconstructedEnvMap
    ))
    {
        requestRecompile();
    }

    if (widget.checkbox(
        "Show SH Grid",
        mbShowAdaptiveGrid
    ))
    {
        requestRecompile();
    }

    if (mbShowAdaptiveGrid)
    {
        if (widget.checkbox(
            "Show Normal Voxels",
            mShowNormalVoxels
        ))
        {
            if (mpProbeVisualizePass)
            {
                mpProbeVisualizePass->setShowNormal(
                    mShowNormalVoxels
                );
            }
        }
        if (widget.checkbox(
            "Show Edge-Added Voxels",
            mbShowEdgeAddedVoxels
        ))
        {
            if (mpProbeVisualizePass)
            {
                mpProbeVisualizePass->setShowEdgeAdded(
                    mbShowEdgeAddedVoxels
                );
            }
        }

        if (widget.checkbox(
            "Draw Leaf Only",
            mbDrawLeafOnly
        ))
        {
            if (mpProbeVisualizePass)
            {
                mpProbeVisualizePass->setDrawLeafOnly(
                    mbDrawLeafOnly
                );
            }
        }

        if (auto g =
            widget.group("Octree Levels", true))
        {
            for (int i = 0; i < 8; ++i)
            {
                std::string label =
                    "Level " + std::to_string(i);

                if (g.checkbox(
                    label.c_str(),
                    mVisLevels[i]
                ))
                {
                    if (mpProbeVisualizePass)
                    {
                        mpProbeVisualizePass->toggleLevel(
                            i,
                            mVisLevels[i]
                        );
                    }
                }
            }
        }
    }
}

void PrecomputeSHCoefficients::setScene(RenderContext* pRenderContext, const ref<Scene>& pScene)
{
    // Set new scene.
    mpScene = pScene;
    if (mpScene)
    {
        
            generateUniformSphereDirSamples(numSamplesPerProbe, samplingDirs); // polar z-up
            //generateProgressiveSphereDirSamples(kMaxSamplesPerProbe, samplingDirs);

            mpProbeDirSamplesBuffer = mpDevice->createStructuredBuffer(
                sizeof(ProbeDirSample), numSamplesPerProbe, ResourceBindFlags::ShaderResource, MemoryType::DeviceLocal, samplingDirs.data()
            );
            mpProbeDirSamplesBuffer->setName("Probe Dir Samples");
            initSHBasisGradientAndHessianTables(samplingDirs);

            //exportIrradianceFieldErrorComparison();
            //exportAnalyticErrorVsDistanceTest();
            //exportResidualTopologyRegionErrorComparison();
           // program
           ProgramDesc desc;
           desc.addShaderModules(mpScene->getShaderModules());
           desc.addShaderLibrary(kShaderFile)
               .vsEntry("vsMain")  // Vertex shader entry point
               .psEntry("psMain"); // Pixel shader entry point;
           mpProgram = Program::create(mpDevice, desc, mpScene->getSceneDefines());
           mpVars = ProgramVars::create(mpDevice, mpProgram->getReflector());

           // rasterizer state
           RasterizerState::Desc rasterDesc;
           rasterDesc.setFillMode(RasterizerState::FillMode::Solid);
           rasterDesc.setCullMode(RasterizerState::CullMode::None);
           rasterDesc.setDepthBias(10000, 1.0f);
           mpRasterState = RasterizerState::create(rasterDesc);

           // default depth stencil state
           DepthStencilState::Desc dsDesc;
           ref<DepthStencilState> pDsState = DepthStencilState::create(dsDesc);

           mpGraphicsState = GraphicsState::create(mpDevice);
           mpGraphicsState->setProgram(mpProgram);
           mpGraphicsState->setRasterizerState(mpRasterState);
           mpGraphicsState->setFbo(mpFbo);
           mpGraphicsState->setDepthStencilState(pDsState);

           //mpFullScreenPass = FullScreenPass::create(mpDevice, kEnvMapShaderFile, mpScene->getSceneDefines(), 0, "vsMain");
           mpProbeVisualizePass = ProbeVisualizePass::create(mpDevice, mpScene->getSceneDefines());

           //Re-apply the visibility masks from our member variables
           for (int i = 0; i < 8; ++i)
           {
               mpProbeVisualizePass->toggleLevel(i, mVisLevels[i]);
           }
           // Restore Leaf Only
           mpProbeVisualizePass->setDrawLeafOnly(mbDrawLeafOnly);
           mpProbeVisualizePass->setShowEdgeAdded(
               mbShowEdgeAddedVoxels
           );
           mpProbeVisualizePass->setShowNormal(
               mShowNormalVoxels
           );

#if CURRENT_PROBE_MODE == PROBE_MODE_ADAPTIVE
           mAdaptiveProbeVolume = AdaptiveProbeVolume::create(mpDevice);
           if (!mNeedRebuildProbeVolume) {
               mAdaptiveProbeVolume->loadFromFile(loadFromFileName);
               mAdaptiveProbeVolume->uploadToGPU();
               mpProbeVisualizePass->setVolumeData(mAdaptiveProbeVolume->getProbes());
           }
#else
           mUniformProbeVolume = UniformProbeVolume::create(mpDevice);
           if (!mNeedRebuildProbeVolume) {
               mUniformProbeVolume->loadFromFile(loadFromFileName);
               mUniformProbeVolume->uploadToGPU();
               mpProbeVisualizePass->setUniformVolumeData(
                   mUniformProbeVolume->getMinPoint(),
                   mUniformProbeVolume->getCellSize(),
                   mUniformProbeVolume->getCellResolution()
               );
           }
#endif
            ProgramDesc rtProgDesc;
            rtProgDesc.addShaderModules(mpScene->getShaderModules());
            rtProgDesc.addShaderLibrary(kProbeSamplingFile);
            rtProgDesc.setMaxTraceRecursionDepth(3); // 1 for calling TraceRay from RayGen, 1 for calling it from the
                                                     // primary-ray ClosestHit shader for reflections, 1 for reflection ray
                                                     // tracing a shadow ray
            rtProgDesc.setMaxPayloadSize(128);        // The largest ray payload struct (PrimaryRayData) is 24 bytes. The payload size
                                                     // should be set as small as possible for maximum performance.
            rtProgDesc.setMaxAttributeSize(8);
            // Add global type conformances.
            rtProgDesc.addTypeConformances(mpScene->getTypeConformances());

            ref<RtBindingTable> sbt = RtBindingTable::create(2, 2, mpScene->getGeometryCount());
            sbt->setRayGen(rtProgDesc.addRayGen("rayGen"));
            sbt->setMiss(0, rtProgDesc.addMiss("primaryMiss"));
             sbt->setMiss(1, rtProgDesc.addMiss("shadowMiss"));
            auto primary = rtProgDesc.addHitGroup("primaryClosestHit");
             auto shadow = rtProgDesc.addHitGroup("", "shadowAnyHit");

            sbt->setHitGroup(0, mpScene->getGeometryIDs(Scene::GeometryType::TriangleMesh), primary);
            sbt->setHitGroup(1, mpScene->getGeometryIDs(Scene::GeometryType::TriangleMesh), shadow);

            const auto& pLights = mpScene->getILightCollection(pRenderContext); //REMARK weird design that light collection is createdupon first call to this.
            if (mpScene->useEmissiveLights())
            {
                if (!mpEmissiveSampler)
                {
                    FALCOR_ASSERT(pLights && pLights->getActiveLightCount(pRenderContext) > 0);
                    FALCOR_ASSERT(!mpEmissiveSampler);

                    switch (mEmissiveSamplerType)
                    {
                        case EmissiveLightSamplerType::Uniform: // use uniform sampling as default for now
                            mpEmissiveSampler =
                                std::make_unique<EmissiveUniformSampler>(pRenderContext, mpScene->getILightCollection(pRenderContext));
                            break;
                        case EmissiveLightSamplerType::LightBVH:
                            mpEmissiveSampler = std::make_unique<LightBVHSampler>(
                                pRenderContext, mpScene->getILightCollection(pRenderContext), mLightBVHOptions
                            );
                            break;
                        case EmissiveLightSamplerType::Power:
                            mpEmissiveSampler =
                                std::make_unique<EmissivePowerSampler>(pRenderContext, mpScene->getILightCollection(pRenderContext));
                            break;
                        default:
                            FALCOR_THROW("Unknown emissive light sampler type");
                    }
                }
            }

           mpRtProgram = Program::create(mpDevice, rtProgDesc, mpScene->getSceneDefines());
        
            if (mpEmissiveSampler)
            {
                auto defines = mpEmissiveSampler->getDefines();
                mpRtProgram->addDefines(defines);
            }

            DefineList lightRelatedDefines;
            lightRelatedDefines.add("USE_ANALYTIC_LIGHTS", mpScene->useAnalyticLights() ? "1" : "0");
            lightRelatedDefines.add("USE_EMISSIVE_LIGHTS", mpScene->useEmissiveLights() ? "1" : "0");

            mpRtProgram->addDefines(lightRelatedDefines);

            mpRtVars = RtProgramVars::create(mpDevice, mpRtProgram, sbt);
            //REMARK :  set all materials to diffuse for SH testing
            auto allMat = pScene->getMaterials();

            for (auto& pMat : allMat)
            {
                // STEP 1: Handle the Base properties (Legacy & Common)
                // Since StandardMaterial inherits BasicMaterial, this runs for EVERYONE.
                auto pBasicMat = pMat->toBasicMaterial();

                if (pBasicMat)
                {
                    // 1. Kill the Specular Color / Shininess
                    // For Legacy OBJ: This makes it matte.
                    // For PBR: This ensures the "F0" (Reflectivity at 0 degrees) is black.
                    pBasicMat->setSpecularParams(float4(0.0f));

                    // 2. Kill Transmission (Glass/Ghosting)
                    pBasicMat->setTransmissionColor(float3(0.0f));
                    pBasicMat->setSpecularTransmission(0.0f);
                    pBasicMat->setDiffuseTransmission(0.0f);
                }

                // STEP 2: Handle the PBR-specific properties
                // This ONLY runs if the material is actually the modern StandardMaterial type.
                StandardMaterial* pStdMat = dynamic_cast<StandardMaterial*>(pMat.get());

                if (pStdMat)
                {
                    // 3. Force PBR Roughness (The most important setting for modern renderers)
                    pStdMat->setRoughness(1.0f);   // 1.0 = Chalk
                    pStdMat->setMetallic(0.0f);    // 0.0 = Dielectric
                    pStdMat->setSpecularTransmission(0.0f);
                    pStdMat->setTransmissionColor(float3(0.0f));
                }
            }
    }
}

// ============================================================================
// Residual topology region / decision quality evaluation
// ============================================================================

struct PairwiseRegionErrorStats
{
    std::vector<float> hessianAPEValues;
    std::vector<float> residualAPEValues;
    std::vector<float> deltaAPEValues;

    double sumHessianAbsErr = 0.0;
    double sumResidualAbsErr = 0.0;
    double sumDeltaAbsErr = 0.0;

    double sumHessianAPE = 0.0;
    double sumResidualAPE = 0.0;
    double sumDeltaAPE = 0.0;

    uint64_t sampleCount = 0;

    uint64_t residualBetterCount = 0;
    uint64_t residualWorseCount = 0;
    uint64_t residualEqualCount = 0;

    double sumPositiveDeltaAbsErr = 0.0;
    double sumNegativeDeltaAbsErr = 0.0;

    double sumPositiveDeltaAPE = 0.0;
    double sumNegativeDeltaAPE = 0.0;
};

struct DecisionRegionErrorStats
{
    // One entry corresponds to one topology-decision region.
    //
    // Prune decision:
    //     key = residual leaf index.
    //
    // Refine decision:
    //     key = original Hessian leaf index.
    //
    // Each decision accumulates all query samples inside that decision region.
    uint64_t sampleCount = 0;

    double sumHessianAbsErr = 0.0;
    double sumResidualAbsErr = 0.0;
    double sumDeltaAbsErr = 0.0;

    double sumHessianAPE = 0.0;
    double sumResidualAPE = 0.0;
    double sumDeltaAPE = 0.0;

    // ------------------------------------------------------------
// Leaf-level diagnostics.
//
// For prune:
//   decision key = residual leaf.
//   residual level is usually fixed.
//   hessian level may vary inside that coarser residual leaf.
//
// For refine:
//   decision key = hessian leaf.
//   hessian level is usually fixed.
//   residual level may vary inside that refined region.
// ------------------------------------------------------------
    bool hasLevelData = false;

    int minHessianLevel = 0;
    int maxHessianLevel = 0;
    int minResidualLevel = 0;
    int maxResidualLevel = 0;

    double sumHessianLevel = 0.0;
    double sumResidualLevel = 0.0;
    double sumLevelDeltaResidualMinusHessian = 0.0;
};

enum class DecisionQualityBucket
{
    ResidualBetterGt20,
    ResidualBetter10To20,
    ResidualBetter5To10,
    Within5Percent,
    ResidualWorse5To10,
    ResidualWorse10To20,
    ResidualWorseGt20
};

struct DecisionQualityBucketCounts
{
    uint64_t betterGt20 = 0;
    uint64_t better10To20 = 0;
    uint64_t better5To10 = 0;
    uint64_t within5 = 0;
    uint64_t worse5To10 = 0;
    uint64_t worse10To20 = 0;
    uint64_t worseGt20 = 0;

    void add(DecisionQualityBucket bucket)
    {
        switch (bucket)
        {
        case DecisionQualityBucket::ResidualBetterGt20:
            betterGt20++;
            break;

        case DecisionQualityBucket::ResidualBetter10To20:
            better10To20++;
            break;

        case DecisionQualityBucket::ResidualBetter5To10:
            better5To10++;
            break;

        case DecisionQualityBucket::Within5Percent:
            within5++;
            break;

        case DecisionQualityBucket::ResidualWorse5To10:
            worse5To10++;
            break;

        case DecisionQualityBucket::ResidualWorse10To20:
            worse10To20++;
            break;

        case DecisionQualityBucket::ResidualWorseGt20:
            worseGt20++;
            break;
        }
    }
};

static double meanFloatVectorForDecisionStats(
    const std::vector<float>& values
)
{
    if (values.empty())
        return 0.0;

    double sum = 0.0;

    for (float v : values)
        sum += double(v);

    return sum / double(values.size());
}

static double computeDecisionDeltaAbsErrPercent(
    double meanHessianAbsErr,
    double meanResidualAbsErr,
    float absTol
)
{
    // Positive:
    //     residual has lower error than Hessian.
    //
    // Negative:
    //     residual has higher error than Hessian.
    //
    // The denominator is clamped by absTol so that very small Hessian errors
    // do not produce unstable percentages.
    double denom =
        std::max(std::abs(meanHessianAbsErr), double(absTol));

    return
        100.0 *
        (meanHessianAbsErr - meanResidualAbsErr) /
        denom;
}

static DecisionQualityBucket classifyDecisionBucket(
    double deltaAbsErrPercent
)
{
    if (deltaAbsErrPercent >= 20.0)
        return DecisionQualityBucket::ResidualBetterGt20;

    if (deltaAbsErrPercent >= 10.0)
        return DecisionQualityBucket::ResidualBetter10To20;

    if (deltaAbsErrPercent >= 5.0)
        return DecisionQualityBucket::ResidualBetter5To10;

    if (deltaAbsErrPercent >= -5.0)
        return DecisionQualityBucket::Within5Percent;

    if (deltaAbsErrPercent >= -10.0)
        return DecisionQualityBucket::ResidualWorse5To10;

    if (deltaAbsErrPercent >= -20.0)
        return DecisionQualityBucket::ResidualWorse10To20;

    return DecisionQualityBucket::ResidualWorseGt20;
}

static std::string bucketNameForPrune(
    DecisionQualityBucket bucket
)
{
    switch (bucket)
    {
    case DecisionQualityBucket::ResidualBetterGt20:
        return "prune_beneficial_gt20pct";

    case DecisionQualityBucket::ResidualBetter10To20:
        return "prune_beneficial_10_20pct";

    case DecisionQualityBucket::ResidualBetter5To10:
        return "prune_beneficial_5_10pct";

    case DecisionQualityBucket::Within5Percent:
        return "prune_within_5pct_tolerance";

    case DecisionQualityBucket::ResidualWorse5To10:
        return "prune_small_loss_5_10pct";

    case DecisionQualityBucket::ResidualWorse10To20:
        return "prune_higher_loss_10_20pct";

    case DecisionQualityBucket::ResidualWorseGt20:
        return "prune_severe_loss_gt20pct";
    }

    return "unknown_prune_bucket";
}

static std::string bucketNameForRefine(
    DecisionQualityBucket bucket
)
{
    switch (bucket)
    {
    case DecisionQualityBucket::ResidualBetterGt20:
        return "refine_strong_improvement_gt20pct";

    case DecisionQualityBucket::ResidualBetter10To20:
        return "refine_useful_improvement_10_20pct";

    case DecisionQualityBucket::ResidualBetter5To10:
        return "refine_small_improvement_5_10pct";

    case DecisionQualityBucket::Within5Percent:
        return "refine_neutral_within_5pct";

    case DecisionQualityBucket::ResidualWorse5To10:
        return "refine_small_loss_5_10pct";

    case DecisionQualityBucket::ResidualWorse10To20:
        return "refine_harmful_loss_10_20pct";

    case DecisionQualityBucket::ResidualWorseGt20:
        return "refine_severe_loss_gt20pct";
    }

    return "unknown_refine_bucket";
}

static std::string decisionClassName(
    const std::string& decisionType,
    DecisionQualityBucket bucket
)
{
    if (decisionType == "prune")
        return bucketNameForPrune(bucket);

    return bucketNameForRefine(bucket);
}

static void addPairwiseRegionErrorSample(
    PairwiseRegionErrorStats& stats,
    float3 hessianIrr,
    float3 residualIrr,
    float3 refIrr,
    float relEpsilon
)
{
    float yRef = luminance709(refIrr);
    float yHessian = luminance709(hessianIrr);
    float yResidual = luminance709(residualIrr);

    float absErrHessian =
        std::abs(yHessian - yRef);

    float absErrResidual =
        std::abs(yResidual - yRef);

    float denom =
        std::max(std::abs(yRef), relEpsilon);

    float apeHessian =
        absErrHessian / denom * 100.0f;

    float apeResidual =
        absErrResidual / denom * 100.0f;

    if (!std::isfinite(apeHessian) ||
        !std::isfinite(apeResidual))
    {
        return;
    }

    float deltaAPE =
        apeHessian - apeResidual;

    float deltaAbsErr =
        absErrHessian - absErrResidual;

    stats.hessianAPEValues.push_back(apeHessian);
    stats.residualAPEValues.push_back(apeResidual);
    stats.deltaAPEValues.push_back(deltaAPE);

    stats.sumHessianAPE += double(apeHessian);
    stats.sumResidualAPE += double(apeResidual);
    stats.sumDeltaAPE += double(deltaAPE);

    stats.sumHessianAbsErr += double(absErrHessian);
    stats.sumResidualAbsErr += double(absErrResidual);
    stats.sumDeltaAbsErr += double(deltaAbsErr);

    stats.sampleCount++;

    const float equalEps = 1e-8f;

    if (absErrResidual + equalEps < absErrHessian)
    {
        stats.residualBetterCount++;

        stats.sumPositiveDeltaAbsErr += double(deltaAbsErr);
        stats.sumPositiveDeltaAPE += double(deltaAPE);
    }
    else if (absErrResidual > absErrHessian + equalEps)
    {
        stats.residualWorseCount++;

        stats.sumNegativeDeltaAbsErr += double(deltaAbsErr);
        stats.sumNegativeDeltaAPE += double(deltaAPE);
    }
    else
    {
        stats.residualEqualCount++;
    }
}

static void addDecisionRegionErrorSample(
    DecisionRegionErrorStats& stats,
    float3 hessianIrr,
    float3 residualIrr,
    float3 refIrr,
    float relEpsilon,
    int hessianLevel,
    int residualLevel
)
{
    float yRef = luminance709(refIrr);
    float yHessian = luminance709(hessianIrr);
    float yResidual = luminance709(residualIrr);

    float absErrHessian =
        std::abs(yHessian - yRef);

    float absErrResidual =
        std::abs(yResidual - yRef);

    float denom =
        std::max(std::abs(yRef), relEpsilon);

    float apeHessian =
        absErrHessian / denom * 100.0f;

    float apeResidual =
        absErrResidual / denom * 100.0f;

    if (!std::isfinite(apeHessian) ||
        !std::isfinite(apeResidual))
    {
        return;
    }

    float deltaAPE =
        apeHessian - apeResidual;

    float deltaAbsErr =
        absErrHessian - absErrResidual;

    stats.sampleCount++;
    // ------------------------------------------------------------
// Leaf-level statistics.
// These are sample-weighted. Because each position is evaluated
// with multiple normals, this is also normal-sample weighted.
// ------------------------------------------------------------
    if (!stats.hasLevelData)
    {
        stats.hasLevelData = true;

        stats.minHessianLevel = hessianLevel;
        stats.maxHessianLevel = hessianLevel;

        stats.minResidualLevel = residualLevel;
        stats.maxResidualLevel = residualLevel;
    }
    else
    {
        stats.minHessianLevel =
            std::min(stats.minHessianLevel, hessianLevel);

        stats.maxHessianLevel =
            std::max(stats.maxHessianLevel, hessianLevel);

        stats.minResidualLevel =
            std::min(stats.minResidualLevel, residualLevel);

        stats.maxResidualLevel =
            std::max(stats.maxResidualLevel, residualLevel);
    }

    stats.sumHessianLevel +=
        double(hessianLevel);

    stats.sumResidualLevel +=
        double(residualLevel);

    stats.sumLevelDeltaResidualMinusHessian +=
        double(residualLevel - hessianLevel);

    stats.sumHessianAbsErr += double(absErrHessian);
    stats.sumResidualAbsErr += double(absErrResidual);
    stats.sumDeltaAbsErr += double(deltaAbsErr);

    stats.sumHessianAPE += double(apeHessian);
    stats.sumResidualAPE += double(apeResidual);
    stats.sumDeltaAPE += double(deltaAPE);
}

static void writePairwiseRegionErrorRow(
    std::ofstream& csv,
    const std::string& region,
    const PairwiseRegionErrorStats& stats,
    uint64_t totalSampleCount
)
{
    double sampleSharePercent =
        totalSampleCount > 0
        ? 100.0 * double(stats.sampleCount) / double(totalSampleCount)
        : 0.0;

    double meanHessianMAPE =
        stats.sampleCount > 0
        ? stats.sumHessianAPE / double(stats.sampleCount)
        : 0.0;

    double meanResidualMAPE =
        stats.sampleCount > 0
        ? stats.sumResidualAPE / double(stats.sampleCount)
        : 0.0;

    double meanDeltaMAPE =
        stats.sampleCount > 0
        ? stats.sumDeltaAPE / double(stats.sampleCount)
        : 0.0;

    double medianDeltaMAPE =
        percentileVector(stats.deltaAPEValues, 50.0);

    double p05DeltaMAPE =
        percentileVector(stats.deltaAPEValues, 5.0);

    double p95DeltaMAPE =
        percentileVector(stats.deltaAPEValues, 95.0);

    double meanHessianAbsErr =
        stats.sampleCount > 0
        ? stats.sumHessianAbsErr / double(stats.sampleCount)
        : 0.0;

    double meanResidualAbsErr =
        stats.sampleCount > 0
        ? stats.sumResidualAbsErr / double(stats.sampleCount)
        : 0.0;

    double meanDeltaAbsErr =
        stats.sampleCount > 0
        ? stats.sumDeltaAbsErr / double(stats.sampleCount)
        : 0.0;

    double residualBetterRate =
        stats.sampleCount > 0
        ? 100.0 * double(stats.residualBetterCount) / double(stats.sampleCount)
        : 0.0;

    double residualWorseRate =
        stats.sampleCount > 0
        ? 100.0 * double(stats.residualWorseCount) / double(stats.sampleCount)
        : 0.0;

    double residualEqualRate =
        stats.sampleCount > 0
        ? 100.0 * double(stats.residualEqualCount) / double(stats.sampleCount)
        : 0.0;

    double meanDeltaAPE_WhenResidualBetter =
        stats.residualBetterCount > 0
        ? stats.sumPositiveDeltaAPE / double(stats.residualBetterCount)
        : 0.0;

    double meanDeltaAPE_WhenResidualWorse =
        stats.residualWorseCount > 0
        ? stats.sumNegativeDeltaAPE / double(stats.residualWorseCount)
        : 0.0;

    double meanDeltaAbsErr_WhenResidualBetter =
        stats.residualBetterCount > 0
        ? stats.sumPositiveDeltaAbsErr / double(stats.residualBetterCount)
        : 0.0;

    double meanDeltaAbsErr_WhenResidualWorse =
        stats.residualWorseCount > 0
        ? stats.sumNegativeDeltaAbsErr / double(stats.residualWorseCount)
        : 0.0;

    csv
        << region << ","
        << stats.sampleCount << ","
        << sampleSharePercent << ","

        << meanHessianMAPE << ","
        << meanResidualMAPE << ","
        << meanDeltaMAPE << ","
        << medianDeltaMAPE << ","
        << p05DeltaMAPE << ","
        << p95DeltaMAPE << ","

        << meanHessianAbsErr << ","
        << meanResidualAbsErr << ","
        << meanDeltaAbsErr << ","

        << residualBetterRate << ","
        << residualWorseRate << ","
        << residualEqualRate << ","

        << meanDeltaAPE_WhenResidualBetter << ","
        << meanDeltaAPE_WhenResidualWorse << ","
        << meanDeltaAbsErr_WhenResidualBetter << ","
        << meanDeltaAbsErr_WhenResidualWorse
        << "\n";
}

static void writeDecisionDetailsRows(
    std::ofstream& csv,
    const std::string& decisionType,
    const std::unordered_map<int, DecisionRegionErrorStats>& decisionStats,
    float absTol
)
{
    for (const auto& kv : decisionStats)
    {
        int decisionKey = kv.first;
        const DecisionRegionErrorStats& s = kv.second;

        if (s.sampleCount == 0)
            continue;

        double meanHessianAbsErr =
            s.sumHessianAbsErr / double(s.sampleCount);

        double meanResidualAbsErr =
            s.sumResidualAbsErr / double(s.sampleCount);

        double meanDeltaAbsErr =
            s.sumDeltaAbsErr / double(s.sampleCount);

        double meanHessianMAPE =
            s.sumHessianAPE / double(s.sampleCount);

        double meanResidualMAPE =
            s.sumResidualAPE / double(s.sampleCount);

        double meanDeltaMAPE =
            s.sumDeltaAPE / double(s.sampleCount);

        double deltaAbsErrPercent =
            computeDecisionDeltaAbsErrPercent(
                meanHessianAbsErr,
                meanResidualAbsErr,
                absTol
            );

        DecisionQualityBucket bucket =
            classifyDecisionBucket(deltaAbsErrPercent);

        std::string decisionClass =
            decisionClassName(decisionType, bucket);
        double meanHessianLevel =
            s.sampleCount > 0
            ? s.sumHessianLevel / double(s.sampleCount)
            : 0.0;

        double meanResidualLevel =
            s.sampleCount > 0
            ? s.sumResidualLevel / double(s.sampleCount)
            : 0.0;

        double meanLevelDelta =
            s.sampleCount > 0
            ? s.sumLevelDeltaResidualMinusHessian / double(s.sampleCount)
            : 0.0;
        csv
            << decisionType << ","
            << decisionKey << ","
            << decisionClass << ","
            << s.sampleCount << ","

            << s.minHessianLevel << ","
            << s.maxHessianLevel << ","
            << s.minResidualLevel << ","
            << s.maxResidualLevel << ","
            << meanHessianLevel << ","
            << meanResidualLevel << ","
            << meanLevelDelta << ","

            << meanHessianAbsErr << ","
            << meanResidualAbsErr << ","
            << meanDeltaAbsErr << ","
            << deltaAbsErrPercent << ","

            << meanHessianMAPE << ","
            << meanResidualMAPE << ","
            << meanDeltaMAPE
            << "\n";
    }
}

static void writeDecisionQualitySummaryRow(
    std::ofstream& csv,
    const std::string& decisionType,
    const std::unordered_map<int, DecisionRegionErrorStats>& decisionStats,
    uint64_t totalNormalSamples,
    float absTol
)
{
    uint64_t decisionCount = 0;
    uint64_t totalDecisionSamples = 0;

    DecisionQualityBucketCounts counts;

    std::vector<float> decisionDeltaAbsErrValues;
    std::vector<float> decisionDeltaAPEValues;
    std::vector<float> decisionDeltaAbsErrPercentValues;

    double sumDecisionHessianAbsErr = 0.0;
    double sumDecisionResidualAbsErr = 0.0;
    double sumDecisionDeltaAbsErr = 0.0;

    double sumDecisionHessianMAPE = 0.0;
    double sumDecisionResidualMAPE = 0.0;
    double sumDecisionDeltaMAPE = 0.0;

    for (const auto& kv : decisionStats)
    {
        const DecisionRegionErrorStats& s = kv.second;

        if (s.sampleCount == 0)
            continue;

        double meanHessianAbsErr =
            s.sumHessianAbsErr / double(s.sampleCount);

        double meanResidualAbsErr =
            s.sumResidualAbsErr / double(s.sampleCount);

        double meanDeltaAbsErr =
            s.sumDeltaAbsErr / double(s.sampleCount);

        double meanHessianMAPE =
            s.sumHessianAPE / double(s.sampleCount);

        double meanResidualMAPE =
            s.sumResidualAPE / double(s.sampleCount);

        double meanDeltaMAPE =
            s.sumDeltaAPE / double(s.sampleCount);

        double deltaAbsErrPercent =
            computeDecisionDeltaAbsErrPercent(
                meanHessianAbsErr,
                meanResidualAbsErr,
                absTol
            );

        DecisionQualityBucket bucket =
            classifyDecisionBucket(deltaAbsErrPercent);

        counts.add(bucket);

        decisionCount++;
        totalDecisionSamples += s.sampleCount;

        sumDecisionHessianAbsErr += meanHessianAbsErr;
        sumDecisionResidualAbsErr += meanResidualAbsErr;
        sumDecisionDeltaAbsErr += meanDeltaAbsErr;

        sumDecisionHessianMAPE += meanHessianMAPE;
        sumDecisionResidualMAPE += meanResidualMAPE;
        sumDecisionDeltaMAPE += meanDeltaMAPE;

        decisionDeltaAbsErrValues.push_back(float(meanDeltaAbsErr));
        decisionDeltaAPEValues.push_back(float(meanDeltaMAPE));
        decisionDeltaAbsErrPercentValues.push_back(float(deltaAbsErrPercent));
    }

    double sampleSharePercent =
        totalNormalSamples > 0
        ? 100.0 * double(totalDecisionSamples) / double(totalNormalSamples)
        : 0.0;

    double meanSamplesPerDecision =
        decisionCount > 0
        ? double(totalDecisionSamples) / double(decisionCount)
        : 0.0;

    auto percentOfDecisions =
        [decisionCount](uint64_t v) -> double
        {
            return decisionCount > 0
                ? 100.0 * double(v) / double(decisionCount)
                : 0.0;
        };

    double meanDecisionHessianAbsErr =
        decisionCount > 0
        ? sumDecisionHessianAbsErr / double(decisionCount)
        : 0.0;

    double meanDecisionResidualAbsErr =
        decisionCount > 0
        ? sumDecisionResidualAbsErr / double(decisionCount)
        : 0.0;

    double meanDecisionDeltaAbsErr =
        decisionCount > 0
        ? sumDecisionDeltaAbsErr / double(decisionCount)
        : 0.0;

    double meanDecisionHessianMAPE =
        decisionCount > 0
        ? sumDecisionHessianMAPE / double(decisionCount)
        : 0.0;

    double meanDecisionResidualMAPE =
        decisionCount > 0
        ? sumDecisionResidualMAPE / double(decisionCount)
        : 0.0;

    double meanDecisionDeltaMAPE =
        decisionCount > 0
        ? sumDecisionDeltaMAPE / double(decisionCount)
        : 0.0;

    uint64_t residualBetterDecisionCount =
        counts.betterGt20 +
        counts.better10To20 +
        counts.better5To10;

    uint64_t withinToleranceDecisionCount =
        counts.within5;

    uint64_t residualWorseDecisionCount =
        counts.worse5To10 +
        counts.worse10To20 +
        counts.worseGt20;

    // Interpretation helper:
    //
    // Prune:
    //     favorable = residual better or within ±5%.
    //     because prune saves memory.
    //
    // Refine:
    //     favorable = residual better by at least 5%.
    //     because refine spends memory.
    //
    // Raw bucket counts are also written, so the interpretation can be changed later.
    uint64_t favorableDecisionCount = 0;
    uint64_t neutralDecisionCount = 0;
    uint64_t unfavorableDecisionCount = 0;

    if (decisionType == "prune")
    {
        favorableDecisionCount =
            residualBetterDecisionCount +
            withinToleranceDecisionCount;

        neutralDecisionCount = 0;

        unfavorableDecisionCount =
            residualWorseDecisionCount;
    }
    else
    {
        favorableDecisionCount =
            residualBetterDecisionCount;

        neutralDecisionCount =
            withinToleranceDecisionCount;

        unfavorableDecisionCount =
            residualWorseDecisionCount;
    }

    csv
        << decisionType << ","
        << decisionCount << ","
        << totalDecisionSamples << ","
        << sampleSharePercent << ","
        << meanSamplesPerDecision << ","

        << favorableDecisionCount << ","
        << neutralDecisionCount << ","
        << unfavorableDecisionCount << ","

        << percentOfDecisions(favorableDecisionCount) << ","
        << percentOfDecisions(neutralDecisionCount) << ","
        << percentOfDecisions(unfavorableDecisionCount) << ","

        << counts.betterGt20 << ","
        << counts.better10To20 << ","
        << counts.better5To10 << ","
        << counts.within5 << ","
        << counts.worse5To10 << ","
        << counts.worse10To20 << ","
        << counts.worseGt20 << ","

        << percentOfDecisions(counts.betterGt20) << ","
        << percentOfDecisions(counts.better10To20) << ","
        << percentOfDecisions(counts.better5To10) << ","
        << percentOfDecisions(counts.within5) << ","
        << percentOfDecisions(counts.worse5To10) << ","
        << percentOfDecisions(counts.worse10To20) << ","
        << percentOfDecisions(counts.worseGt20) << ","

        << meanDecisionHessianMAPE << ","
        << meanDecisionResidualMAPE << ","
        << meanDecisionDeltaMAPE << ","
        << percentileVector(decisionDeltaAPEValues, 50.0) << ","
        << percentileVector(decisionDeltaAPEValues, 5.0) << ","
        << percentileVector(decisionDeltaAPEValues, 95.0) << ","

        << meanDecisionHessianAbsErr << ","
        << meanDecisionResidualAbsErr << ","
        << meanDecisionDeltaAbsErr << ","
        << percentileVector(decisionDeltaAbsErrValues, 50.0) << ","
        << percentileVector(decisionDeltaAbsErrValues, 5.0) << ","
        << percentileVector(decisionDeltaAbsErrValues, 95.0) << ","

        << meanFloatVectorForDecisionStats(decisionDeltaAbsErrPercentValues) << ","
        << percentileVector(decisionDeltaAbsErrPercentValues, 50.0) << ","
        << percentileVector(decisionDeltaAbsErrPercentValues, 5.0) << ","
        << percentileVector(decisionDeltaAbsErrPercentValues, 95.0)
        << "\n";
}

static double safePercentU64(
    uint64_t count,
    uint64_t total
)
{
    return total > 0
        ? 100.0 * double(count) / double(total)
        : 0.0;
}

static void writeRegionTextSummary(
    std::ofstream& out,
    const std::string& regionName,
    const PairwiseRegionErrorStats& stats,
    uint64_t totalNormalSamples
)
{
    double sampleSharePercent =
        totalNormalSamples > 0
        ? 100.0 * double(stats.sampleCount) / double(totalNormalSamples)
        : 0.0;

    double hessianMAPE =
        stats.sampleCount > 0
        ? stats.sumHessianAPE / double(stats.sampleCount)
        : 0.0;

    double residualMAPE =
        stats.sampleCount > 0
        ? stats.sumResidualAPE / double(stats.sampleCount)
        : 0.0;

    double deltaMAPE =
        stats.sampleCount > 0
        ? stats.sumDeltaAPE / double(stats.sampleCount)
        : 0.0;

    double hessianAbsErr =
        stats.sampleCount > 0
        ? stats.sumHessianAbsErr / double(stats.sampleCount)
        : 0.0;

    double residualAbsErr =
        stats.sampleCount > 0
        ? stats.sumResidualAbsErr / double(stats.sampleCount)
        : 0.0;

    double deltaAbsErr =
        stats.sampleCount > 0
        ? stats.sumDeltaAbsErr / double(stats.sampleCount)
        : 0.0;

    double betterRate =
        safePercentU64(stats.residualBetterCount, stats.sampleCount);

    double worseRate =
        safePercentU64(stats.residualWorseCount, stats.sampleCount);

    double equalRate =
        safePercentU64(stats.residualEqualCount, stats.sampleCount);

    out << "------------------------------------------------------------\n";
    out << " Region: " << regionName << "\n";
    out << "------------------------------------------------------------\n";

    out << " Sample count:                    "
        << stats.sampleCount << "\n";

    out << " Sample share:                    "
        << sampleSharePercent << "%\n";

    out << " Hessian MAPE:                    "
        << hessianMAPE << "\n";

    out << " Residual MAPE:                   "
        << residualMAPE << "\n";

    out << " Delta MAPE (Hessian - Residual): "
        << deltaMAPE << "\n";

    out << " Median delta MAPE:               "
        << percentileVector(stats.deltaAPEValues, 50.0) << "\n";

    out << " P05 delta MAPE:                  "
        << percentileVector(stats.deltaAPEValues, 5.0) << "\n";

    out << " P95 delta MAPE:                  "
        << percentileVector(stats.deltaAPEValues, 95.0) << "\n";

    out << " Hessian mean abs error:          "
        << hessianAbsErr << "\n";

    out << " Residual mean abs error:         "
        << residualAbsErr << "\n";

    out << " Delta mean abs error:            "
        << deltaAbsErr << "\n";

    out << " Residual better samples:         "
        << stats.residualBetterCount
        << " (" << betterRate << "%)\n";

    out << " Residual worse samples:          "
        << stats.residualWorseCount
        << " (" << worseRate << "%)\n";

    out << " Residual equal samples:          "
        << stats.residualEqualCount
        << " (" << equalRate << "%)\n";

    out << "\n";
}

static void writeDecisionQualityTextSummary(
    std::ofstream& out,
    const std::string& decisionType,
    const std::unordered_map<int, DecisionRegionErrorStats>& decisionStats,
    uint64_t totalNormalSamples,
    float absTol
)
{
    uint64_t decisionCount = 0;
    uint64_t totalDecisionSamples = 0;

    DecisionQualityBucketCounts counts;

    std::vector<float> decisionDeltaAbsErrValues;
    std::vector<float> decisionDeltaAPEValues;
    std::vector<float> decisionDeltaAbsErrPercentValues;

    double sumDecisionHessianAbsErr = 0.0;
    double sumDecisionResidualAbsErr = 0.0;
    double sumDecisionDeltaAbsErr = 0.0;

    double sumDecisionHessianMAPE = 0.0;
    double sumDecisionResidualMAPE = 0.0;
    double sumDecisionDeltaMAPE = 0.0;

    for (const auto& kv : decisionStats)
    {
        const DecisionRegionErrorStats& s = kv.second;

        if (s.sampleCount == 0)
            continue;

        double meanHessianAbsErr =
            s.sumHessianAbsErr / double(s.sampleCount);

        double meanResidualAbsErr =
            s.sumResidualAbsErr / double(s.sampleCount);

        double meanDeltaAbsErr =
            s.sumDeltaAbsErr / double(s.sampleCount);

        double meanHessianMAPE =
            s.sumHessianAPE / double(s.sampleCount);

        double meanResidualMAPE =
            s.sumResidualAPE / double(s.sampleCount);

        double meanDeltaMAPE =
            s.sumDeltaAPE / double(s.sampleCount);

        double deltaAbsErrPercent =
            computeDecisionDeltaAbsErrPercent(
                meanHessianAbsErr,
                meanResidualAbsErr,
                absTol
            );

        DecisionQualityBucket bucket =
            classifyDecisionBucket(deltaAbsErrPercent);

        counts.add(bucket);

        decisionCount++;
        totalDecisionSamples += s.sampleCount;

        sumDecisionHessianAbsErr += meanHessianAbsErr;
        sumDecisionResidualAbsErr += meanResidualAbsErr;
        sumDecisionDeltaAbsErr += meanDeltaAbsErr;

        sumDecisionHessianMAPE += meanHessianMAPE;
        sumDecisionResidualMAPE += meanResidualMAPE;
        sumDecisionDeltaMAPE += meanDeltaMAPE;

        decisionDeltaAbsErrValues.push_back(float(meanDeltaAbsErr));
        decisionDeltaAPEValues.push_back(float(meanDeltaMAPE));
        decisionDeltaAbsErrPercentValues.push_back(float(deltaAbsErrPercent));
    }

    uint64_t residualBetterDecisionCount =
        counts.betterGt20 +
        counts.better10To20 +
        counts.better5To10;

    uint64_t withinToleranceDecisionCount =
        counts.within5;

    uint64_t residualWorseDecisionCount =
        counts.worse5To10 +
        counts.worse10To20 +
        counts.worseGt20;

    uint64_t favorableDecisionCount = 0;
    uint64_t neutralDecisionCount = 0;
    uint64_t unfavorableDecisionCount = 0;

    if (decisionType == "prune")
    {
        favorableDecisionCount =
            residualBetterDecisionCount +
            withinToleranceDecisionCount;

        neutralDecisionCount = 0;

        unfavorableDecisionCount =
            residualWorseDecisionCount;
    }
    else
    {
        favorableDecisionCount =
            residualBetterDecisionCount;

        neutralDecisionCount =
            withinToleranceDecisionCount;

        unfavorableDecisionCount =
            residualWorseDecisionCount;
    }

    double sampleSharePercent =
        totalNormalSamples > 0
        ? 100.0 * double(totalDecisionSamples) / double(totalNormalSamples)
        : 0.0;

    double meanSamplesPerDecision =
        decisionCount > 0
        ? double(totalDecisionSamples) / double(decisionCount)
        : 0.0;

    double meanDecisionHessianAbsErr =
        decisionCount > 0
        ? sumDecisionHessianAbsErr / double(decisionCount)
        : 0.0;

    double meanDecisionResidualAbsErr =
        decisionCount > 0
        ? sumDecisionResidualAbsErr / double(decisionCount)
        : 0.0;

    double meanDecisionDeltaAbsErr =
        decisionCount > 0
        ? sumDecisionDeltaAbsErr / double(decisionCount)
        : 0.0;

    double meanDecisionHessianMAPE =
        decisionCount > 0
        ? sumDecisionHessianMAPE / double(decisionCount)
        : 0.0;

    double meanDecisionResidualMAPE =
        decisionCount > 0
        ? sumDecisionResidualMAPE / double(decisionCount)
        : 0.0;

    double meanDecisionDeltaMAPE =
        decisionCount > 0
        ? sumDecisionDeltaMAPE / double(decisionCount)
        : 0.0;

    out << "============================================================\n";
    out << " DECISION QUALITY: " << decisionType << "\n";
    out << "============================================================\n";

    out << " Decision count:                  "
        << decisionCount << "\n";

    out << " Total decision samples:           "
        << totalDecisionSamples << "\n";

    out << " Sample share:                     "
        << sampleSharePercent << "%\n";

    out << " Mean samples per decision:        "
        << meanSamplesPerDecision << "\n";

    out << "\n";

    out << " Favorable decisions:              "
        << favorableDecisionCount
        << " (" << safePercentU64(favorableDecisionCount, decisionCount) << "%)\n";

    out << " Neutral decisions:                "
        << neutralDecisionCount
        << " (" << safePercentU64(neutralDecisionCount, decisionCount) << "%)\n";

    out << " Unfavorable decisions:            "
        << unfavorableDecisionCount
        << " (" << safePercentU64(unfavorableDecisionCount, decisionCount) << "%)\n";

    out << "\n";

    out << " Bucket counts:\n";

    out << "   Residual better >20%:           "
        << counts.betterGt20
        << " (" << safePercentU64(counts.betterGt20, decisionCount) << "%)\n";

    out << "   Residual better 10-20%:         "
        << counts.better10To20
        << " (" << safePercentU64(counts.better10To20, decisionCount) << "%)\n";

    out << "   Residual better 5-10%:          "
        << counts.better5To10
        << " (" << safePercentU64(counts.better5To10, decisionCount) << "%)\n";

    out << "   Within +/-5%:                   "
        << counts.within5
        << " (" << safePercentU64(counts.within5, decisionCount) << "%)\n";

    out << "   Residual worse 5-10%:           "
        << counts.worse5To10
        << " (" << safePercentU64(counts.worse5To10, decisionCount) << "%)\n";

    out << "   Residual worse 10-20%:          "
        << counts.worse10To20
        << " (" << safePercentU64(counts.worse10To20, decisionCount) << "%)\n";

    out << "   Residual worse >20%:            "
        << counts.worseGt20
        << " (" << safePercentU64(counts.worseGt20, decisionCount) << "%)\n";

    out << "\n";

    out << " Decision-weighted error:\n";

    out << "   Mean Hessian MAPE:              "
        << meanDecisionHessianMAPE << "\n";

    out << "   Mean Residual MAPE:             "
        << meanDecisionResidualMAPE << "\n";

    out << "   Mean Delta MAPE:                "
        << meanDecisionDeltaMAPE << "\n";

    out << "   Median Delta MAPE:              "
        << percentileVector(decisionDeltaAPEValues, 50.0) << "\n";

    out << "   P05 Delta MAPE:                 "
        << percentileVector(decisionDeltaAPEValues, 5.0) << "\n";

    out << "   P95 Delta MAPE:                 "
        << percentileVector(decisionDeltaAPEValues, 95.0) << "\n";

    out << "\n";

    out << "   Mean Hessian abs error:         "
        << meanDecisionHessianAbsErr << "\n";

    out << "   Mean Residual abs error:        "
        << meanDecisionResidualAbsErr << "\n";

    out << "   Mean Delta abs error:           "
        << meanDecisionDeltaAbsErr << "\n";

    out << "   Median Delta abs error:         "
        << percentileVector(decisionDeltaAbsErrValues, 50.0) << "\n";

    out << "   P05 Delta abs error:            "
        << percentileVector(decisionDeltaAbsErrValues, 5.0) << "\n";

    out << "   P95 Delta abs error:            "
        << percentileVector(decisionDeltaAbsErrValues, 95.0) << "\n";

    out << "\n";

    out << "   Mean Delta abs error percent:   "
        << meanFloatVectorForDecisionStats(decisionDeltaAbsErrPercentValues) << "%\n";

    out << "   Median Delta abs error percent: "
        << percentileVector(decisionDeltaAbsErrPercentValues, 50.0) << "%\n";

    out << "   P05 Delta abs error percent:    "
        << percentileVector(decisionDeltaAbsErrPercentValues, 5.0) << "%\n";

    out << "   P95 Delta abs error percent:    "
        << percentileVector(decisionDeltaAbsErrPercentValues, 95.0) << "%\n";

    out << "\n";
}

static bool writeResidualTopologyMainStatsText(
    const std::string& outTxtFile,
    const std::string& referenceFile,
    const std::string& hessianFile,
    const std::string& residualFile,
    uint32_t queryRes,
    uint64_t totalNormalSamples,
    uint64_t failedTopologyLookup,
    const PairwiseRegionErrorStats& allStats,
    const PairwiseRegionErrorStats& sameLevelStats,
    const PairwiseRegionErrorStats& prunedStats,
    const PairwiseRegionErrorStats& addedStats,
    const std::unordered_map<int, DecisionRegionErrorStats>& pruneDecisionStats,
    const std::unordered_map<int, DecisionRegionErrorStats>& refineDecisionStats,
    float decisionAbsTol
)
{
    std::ofstream out(outTxtFile);

    if (!out)
        return false;

    out << std::fixed << std::setprecision(8);

    out << "============================================================\n";
    out << " RESIDUAL TOPOLOGY MAIN STATISTICS\n";
    out << "============================================================\n";

    out << " Reference file:                  " << referenceFile << "\n";
    out << " Hessian file:                    " << hessianFile << "\n";
    out << " Residual file:                   " << residualFile << "\n";
    out << " Query resolution:                " << queryRes << "^3\n";
    out << " Total normal samples:            " << totalNormalSamples << "\n";
    out << " Failed topology lookups:         " << failedTopologyLookup << "\n";
    out << " Decision percent abs tolerance:  " << decisionAbsTol << "\n";

    out << "\n";

    out << "============================================================\n";
    out << " SAMPLE-WEIGHTED REGION QUALITY\n";
    out << "============================================================\n\n";

    writeRegionTextSummary(
        out,
        "all",
        allStats,
        totalNormalSamples
    );

    writeRegionTextSummary(
        out,
        "same_level",
        sameLevelStats,
        totalNormalSamples
    );

    writeRegionTextSummary(
        out,
        "pruned_residual_coarser",
        prunedStats,
        totalNormalSamples
    );

    writeRegionTextSummary(
        out,
        "added_residual_finer",
        addedStats,
        totalNormalSamples
    );

    writeDecisionQualityTextSummary(
        out,
        "prune",
        pruneDecisionStats,
        totalNormalSamples,
        decisionAbsTol
    );

    writeDecisionQualityTextSummary(
        out,
        "refine",
        refineDecisionStats,
        totalNormalSamples,
        decisionAbsTol
    );

    out << "============================================================\n";
    out << " INTERPRETATION NOTES\n";
    out << "============================================================\n";

    out << " Delta = Hessian error - Residual error.\n";
    out << " Positive delta means residual is better.\n";
    out << " Negative delta means residual is worse.\n\n";

    out << " For prune decisions:\n";
    out << "   Favorable = residual better OR within +/-5%.\n";
    out << "   This is because pruning saves memory, so small loss is acceptable.\n\n";

    out << " For refine decisions:\n";
    out << "   Favorable = residual better by at least 5%.\n";
    out << "   Neutral = within +/-5%.\n";
    out << "   This is because refinement spends memory and should improve quality.\n\n";

    out << " Bucket thresholds are based on decision-level mean absolute error change.\n";

    return true;
}

void PrecomputeSHCoefficients::exportResidualTopologyRegionErrorComparison()
{
    // ============================================================
    // Input files
    // ============================================================

    //const std::string referenceFile =
    //    "U64CornellShadowBoundaryScene.txt";

    //const std::string hessianFile =
    //    "DirectAbsErr2HessianMetricCornellThinSlab.txt";

    //const std::string egcFile =
    //    "DirectAbsErr2ResidualScaleV2MetricCornellThinSlab.txt";
    //const std::string summaryCsvFile =
    //    "CornellSceneEGCAddedRegionErrorSummary.csv";

    const std::string referenceFile =  "U64DataScene_4096spp.txt";
    const std::string hessianFile =  "DirectAbsErr2N6HessianMetricDataScene4096spp.txt";
    const std::string egcFile =  "DirectAbsErr2N6EdgeGradientMetricDataScene4096spp.txt";
    const std::string summaryCsvFile = "DataSceneEGCAddedRegionErrorSummary.csv";
    // ============================================================
    // Output files
    // ============================================================



    // ============================================================
    // Evaluation settings
    // ============================================================

    // Fixed number of spatial samples along each axis of every
    // EGC-added leaf node.
    constexpr uint32_t samplesPerAxis = 4;


    const std::vector<float3> normalDirections =
        makeSixAxisNormalsWorld();

    SampleGainLossStats gainLossStats;
    // ============================================================
    // Load volumes
    // ============================================================

    logInfo(
        "Loading local EGC-added-region evaluation volumes."
    );

    logInfo("  Reference = " + referenceFile);
    logInfo("  Hessian   = " + hessianFile);
    logInfo("  EGC       = " + egcFile);

    ref<UniformProbeVolume> reference =
        UniformProbeVolume::create(mpDevice);

    reference->loadFromFile(referenceFile);

    ref<AdaptiveProbeVolume> hessianGrid =
        AdaptiveProbeVolume::create(mpDevice);

    hessianGrid->loadFromFile(hessianFile);

    ref<AdaptiveProbeVolume> egcGrid =
        AdaptiveProbeVolume::create(mpDevice);

    egcGrid->loadFromFile(egcFile);
 
    // ============================================================
    // Evaluate EGC-added leaf nodes
    // ============================================================

    uint64_t egcAddedLeafCount = 0;
    uint64_t failedReferenceQueries = 0;
    uint64_t failedHessianQueries = 0;
    uint64_t failedEGCQueries = 0;

    const auto& egcProbes =
        egcGrid->getProbes();

    for (
        int egcProbeIndex = 0;
        egcProbeIndex < int(egcProbes.size());
        ++egcProbeIndex
        )
    {
        const auto& probe =
            egcProbes[egcProbeIndex];

        // Evaluate only final leaves to avoid sampling overlapping
        // parent and child regions.
        if (!probe.isLeaf)
            continue;

        // addedByEdgeAware is propagated to descendants, so this
        // selects the final EGC-added topology.
        if (!probe.addedByEdgeAware)
            continue;

        egcAddedLeafCount++;

        const float3 nodeMin =
            probe.minPoint;

        const float3 nodeMax =
            probe.maxPoint;

        const float3 nodeExtent =
            nodeMax - nodeMin;

        for (
            uint32_t z = 0;
            z < samplesPerAxis;
            ++z
            )
        {
            for (
                uint32_t y = 0;
                y < samplesPerAxis;
                ++y
                )
            {
                for (
                    uint32_t x = 0;
                    x < samplesPerAxis;
                    ++x
                    )
                {
                    // Cell-centered local coordinates.
                    // For samplesPerAxis = 4:
                    // 0.125, 0.375, 0.625, 0.875.
                    const float3 localUVW =
                        float3(
                            (float(x) + 0.5f) /
                            float(samplesPerAxis),

                            (float(y) + 0.5f) /
                            float(samplesPerAxis),

                            (float(z) + 0.5f) /
                            float(samplesPerAxis)
                        );

                    const float3 positionWorld =
                        nodeMin +
                        localUVW * nodeExtent;

                    // Validate that all three volumes contain the
                    // position before evaluating irradiance.
                    const int hessianLeafIndex =
                        hessianGrid->
                        findLeafProbeIndexCPU(
                            positionWorld
                        );

                    if (hessianLeafIndex < 0)
                    {
                        failedHessianQueries++;
                        continue;
                    }

                    const int egcLeafIndex =
                        egcGrid->
                        findLeafProbeIndexCPU(
                            positionWorld
                        );

                    if (egcLeafIndex < 0)
                    {
                        failedEGCQueries++;
                        continue;
                    }

                    // The query should normally return the same EGC
                    // leaf currently being sampled. A different leaf
                    // can occur only from boundary precision.
                    const float3 referenceMin =
                        reference->getMinPoint();

                    const float3 referenceMax =
                        reference->getMaxPoint();

                    const bool insideReference =
                        positionWorld.x >= referenceMin.x &&
                        positionWorld.y >= referenceMin.y &&
                        positionWorld.z >= referenceMin.z &&
                        positionWorld.x <= referenceMax.x &&
                        positionWorld.y <= referenceMax.y &&
                        positionWorld.z <= referenceMax.z;

                    if (!insideReference)
                    {
                        failedReferenceQueries++;
                        continue;
                    }

                    for (
                        const float3& normalWorld :
                        normalDirections
                        )
                    {
                        const float3 referenceIrradiance =
                            reference->
                            evaluateIrradianceHermiteCPU(
                                positionWorld,
                                normalWorld
                            );

                        const float3 hessianIrradiance =
                            hessianGrid->
                            evaluateIrradianceHermiteCPU(
                                positionWorld,
                                normalWorld
                            );

                        const float3 egcIrradiance =
                            egcGrid->
                            evaluateIrradianceHermiteCPU(
                                positionWorld,
                                normalWorld
                            );

                        // Collect gain for this position-normal sample.
                        const float referenceLuminance =
                            luminance709(referenceIrradiance);

                        const float hessianLuminance =
                            luminance709(hessianIrradiance);

                        const float egcLuminance =
                            luminance709(egcIrradiance);

                        gainLossStats.add(
                            std::abs(
                                hessianLuminance -
                                referenceLuminance
                            ),
                            std::abs(
                                egcLuminance -
                                referenceLuminance
                            )
                        );
                    }
                }
            }
        }
    }

    // ============================================================
// Write compact gain/loss summary
// ============================================================

    std::ofstream summaryCsv(summaryCsvFile);

    if (!summaryCsv)
    {
        logError(
            "Failed to open EGC gain/loss summary CSV: " +
            summaryCsvFile
        );

        return;
    }

    summaryCsv << std::fixed << std::setprecision(10);

    summaryCsv
        << "referenceFile,"
        << "hessianFile,"
        << "egcFile,"
        << "egcAddedLeafCount,"
        << "validSampleCount,"

        << "improvedCount,"
        << "improvedPercent,"
        << "meanGain,"
        << "medianGain,"
        << "totalGain,"

        << "worsenedCount,"
        << "worsenedPercent,"
        << "meanLoss,"
        << "medianLoss,"
        << "totalLoss,"

        << "comparableCount,"
        << "comparablePercent,"

        << "netGain,"
        << "netMeanErrorReduction,"
        << "gainLossRatio"
        << "\n";

    summaryCsv
        << referenceFile << ","
        << hessianFile << ","
        << egcFile << ","
        << egcAddedLeafCount << ","
        << gainLossStats.getValidSampleCount() << ","

        << gainLossStats.improvedCount << ","
        << gainLossStats.percent(
            gainLossStats.improvedCount
        ) << ","
        << gainLossStats.getMeanGain() << ","
        << gainLossStats.getMedianGain() << ","
        << gainLossStats.totalGain << ","

        << gainLossStats.worsenedCount << ","
        << gainLossStats.percent(
            gainLossStats.worsenedCount
        ) << ","
        << gainLossStats.getMeanLoss() << ","
        << gainLossStats.getMedianLoss() << ","
        << gainLossStats.totalLoss << ","

        << gainLossStats.comparableCount << ","
        << gainLossStats.percent(
            gainLossStats.comparableCount
        ) << ","

        << gainLossStats.getNetGain() << ","
        << gainLossStats.getNetMeanErrorReduction() << ","
        << gainLossStats.getGainLossRatio()
        << "\n";

    summaryCsv.close();
}
