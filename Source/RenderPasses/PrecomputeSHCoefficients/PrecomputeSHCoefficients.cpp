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
//const int numSamplesPerProbe = 1;
const int verificationRes = 100;
const float verificationH = 0.001f;
const float verificationY = 0.2f;
const float verificationExtent = 0.25f;

//const float ErrorThreshold = 10.0f;
//const float ErrorThreshold = 5.0f;
const float ErrorThreshold =3.0f;//threshold for Erel
//const float ErrorThreshold =1.5f;//threshold for Erel
//const bool useRelativeError = false;
const bool useRelativeError = true;
//const uint3 unifromGridSize = uint3(16, 16, 16);
//const uint3 unifromGridSize = uint3(32, 32, 32);
//const uint3 unifromGridSize = uint3(8, 8, 8);
const uint3 unifromGridSize = uint3(64, 64, 64);
const std::string loadFromFileName = "IndirectUniformGrid32.txt";

//const std::string saveToFileName = "Seeded16DirectAbsErr3SubwayCorridorNoOpen.txt";
//const std::string saveToFileName = "Seeded16DirectAbsErr5SubwayCorridorNoOpen.txt";
//const std::string saveToFileName = "Seeded16DirectAbsErr10SubwayCorridorNoOpen.txt";


const std::string saveToFileName = "Seeded8DirectAbsErr5SubwayCorridorNoOpen.txt";

//const std::string saveToFileName = "DirectUniformGrid32SubwayCorridorNoOpen.txt";
//const std::string saveToFileName = "DirectUniformGrid64SubwayCorridorNoOpen.txt";

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
} // namespace

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
            // 1. Initialize
            //mAdaptiveProbeVolume->startBuild(mpScene, ErrorThreshold, useRelativeError);

            const uint3 seedResolution = uint3(8,8,8); 
            //const uint3 seedResolution = uint3(16, 16, 16); 
            mAdaptiveProbeVolume->startBuildSeeded(mpScene, seedResolution, ErrorThreshold, useRelativeError);

            const uint32_t kMaxCornersPerDispatch = 8192; // tune this

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
                        std::vector<ProbeSampleData> probeSamplingResults;
                        probeSamplingResults.reserve(numSamplesPerProbe);
                        for (int sampleIdx = 0; sampleIdx < numSamplesPerProbe; ++sampleIdx)
                            probeSamplingResults.push_back(allProbeSamplingData[offset + sampleIdx]);

                        float3 xPolar;
                        xPolar.x = pendingProbePositions[probeIdx].z;
                        xPolar.y = pendingProbePositions[probeIdx].x;
                        xPolar.z = pendingProbePositions[probeIdx].y;

                        calculateSHCoeffsGradientsRGBAndHessiansLum(
                            gradsBatch[probeIdx],
                            hessiansBatch[probeIdx],
                            xPolar,
                            probeSamplingResults,
                            samplingDirs
                        );
                        calculateSHCoeffs(coeffsBatch[probeIdx], probeSamplingResults, numSamplesPerProbe);
                    }

                    mAdaptiveProbeVolume->setCornerDataRange(batchStart, coeffsBatch, gradsBatch, hessiansBatch);
                }

                mAdaptiveProbeVolume->finishBatch();
            }

            auto tEnd = clock::now();
            double ms = std::chrono::duration<double, std::milli>(tEnd - tStart).count();
            mAdaptiveProbeVolume->setBuildTimeMs(ms);
            // 3. TODO: upload to GPU for visualization
            mAdaptiveProbeVolume->uploadToGPU();
            mAdaptiveProbeVolume->printDebugInfo("AdaptiveProbeVolumeTextViz.txt");
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

                        int offset = x * numSamplesPerProbe;
                        std::vector<ProbeSampleData> probeSamplingResults;
                        probeSamplingResults.assign(rowSamplingData.begin() + offset, rowSamplingData.begin() + offset + numSamplesPerProbe);

                        // Convert to Polar for SH Math (Z, X, Y mapping)
                        float3 xPolar = float3(rowPositions[x].z, rowPositions[x].x, rowPositions[x].y);

                        std::vector<float3> coeffs;
                        std::vector<GradSHCoeff> grads;

                        // Perform Physics Calculations
                        calculateSHCoeffsGradients(grads, xPolar, probeSamplingResults, samplingDirs);
                        calculateSHCoeffs(coeffs, probeSamplingResults, numSamplesPerProbe);

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

             if (mbShowAdaptiveGrid)
             {
                mpProbeVisualizePass->setCameraData(
                    mpScene->getCamera()->getViewProjMatrix()
                );
               
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

void PrecomputeSHCoefficients::renderUI(Gui::Widgets& widget)
{
    if (widget.checkbox("Show Reconstructed Env Map", mbShowReconstructedEnvMap))
        requestRecompile();
    if (widget.checkbox("Show SH Grid", mbShowAdaptiveGrid))
        requestRecompile();
    // Level Visibility Controls
    if (mbShowAdaptiveGrid)
    {
        // NEW Checkbox
        if (widget.checkbox("Draw Leaf Only", mbDrawLeafOnly))
        {
            if (mpProbeVisualizePass)
                mpProbeVisualizePass->setDrawLeafOnly(mbDrawLeafOnly);
        }

        if (auto g = widget.group("Octree Levels", true))
        {
            for (int i = 0; i < 8; ++i)
            {
                std::string label = "Level " + std::to_string(i);
                if (g.checkbox(label.c_str(), mVisLevels[i]))
                {
                    if (mpProbeVisualizePass)
                        mpProbeVisualizePass->toggleLevel(i, mVisLevels[i]);
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
            mpProbeDirSamplesBuffer = mpDevice->createStructuredBuffer(
                sizeof(ProbeDirSample), numSamplesPerProbe, ResourceBindFlags::ShaderResource, MemoryType::DeviceLocal, samplingDirs.data()
            );
            mpProbeDirSamplesBuffer->setName("Probe Dir Samples");
            initSHBasisGradientAndHessianTables(samplingDirs);

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

void PrecomputeSHCoefficients::runDistanceVoxelLineExperiment(RenderContext* pRenderContext)
{
    // ------------------------------------------------------------
    // Experiment settings
    // ------------------------------------------------------------
    const float3 lightPos = float3(0.f, 5.0f, 0.05f);

    // Fixed voxel size for the whole test.
    const float3 voxelSize = float3(1.0f, 1.0f, 1.0f);
    const float3 halfSize = voxelSize * 0.5f;

    // Sample voxel centers along corridor centerline.
    const uint32_t numVoxels = 24;
    const float3 lineStart = float3(0.f, 0.5f, 1.0f);
    const float zStep = 1.0f;

    std::vector<VoxelLineVoxel> voxels;
    voxels.reserve(numVoxels);

    std::vector<float3> uniqueCornerPositions;
    uniqueCornerPositions.reserve(numVoxels * 8);

    struct CornerKey
    {
        int x, y, z;
        bool operator==(const CornerKey& other) const
        {
            return x == other.x && y == other.y && z == other.z;
        }
    };

    struct CornerKeyHasher
    {
        std::size_t operator()(const CornerKey& k) const
        {
            return ((std::hash<int>()(k.x) ^ (std::hash<int>()(k.y) << 1)) >> 1) ^
                (std::hash<int>()(k.z) << 1);
        }
    };

    std::unordered_map<CornerKey, int, CornerKeyHasher> cornerLookup;
    const float quantizationScale = 10000.0f;

    auto makeCornerKey = [&](const float3& p) -> CornerKey
        {
            return {
                (int)std::floor(p.x * quantizationScale + 0.5f),
                (int)std::floor(p.y * quantizationScale + 0.5f),
                (int)std::floor(p.z * quantizationScale + 0.5f)
            };
        };

    auto addCorner = [&](const float3& p) -> int
        {
            CornerKey key = makeCornerKey(p);
            auto it = cornerLookup.find(key);
            if (it != cornerLookup.end()) return it->second;

            int idx = (int)uniqueCornerPositions.size();
            uniqueCornerPositions.push_back(p);
            cornerLookup[key] = idx;
            return idx;
        };

    // ------------------------------------------------------------
    // Generate virtual voxels and deduplicated corners
    // ------------------------------------------------------------
    for (uint32_t i = 0; i < numVoxels; ++i)
    {
        VoxelLineVoxel voxel{};
        voxel.center = lineStart + float3(0.f, 0.f, zStep * float(i));

        for (int c = 0; c < 8; ++c)
        {
            float3 offset(
                (c & 4) ? halfSize.x : -halfSize.x,
                (c & 2) ? halfSize.y : -halfSize.y,
                (c & 1) ? halfSize.z : -halfSize.z
            );

            float3 cornerPos = voxel.center + offset;
            voxel.cornerIndices[c] = addCorner(cornerPos);
        }

        voxels.push_back(voxel);
    }

    const uint32_t numCorners = (uint32_t)uniqueCornerPositions.size();
    if (numCorners == 0) return;

    // ------------------------------------------------------------
    // Ray trace all unique corners in one batch
    // ------------------------------------------------------------
    mpProbePosBuffer = mpDevice->createStructuredBuffer(
        sizeof(float3),
        numCorners,
        ResourceBindFlags::ShaderResource,
        MemoryType::DeviceLocal,
        uniqueCornerPositions.data()
    );

    mpProbeSamplingResultBuffer = mpDevice->createStructuredBuffer(
        sizeof(ProbeSampleData),
        numSamplesPerProbe * numCorners,
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

    mpScene->raytrace(pRenderContext, mpRtProgram.get(), mpRtVars, uint3(numSamplesPerProbe, numCorners, 1));

    std::vector<ProbeSampleData> allSamplingData(numSamplesPerProbe * numCorners);
    mpProbeSamplingResultBuffer->getBlob(
        allSamplingData.data(),
        0,
        numSamplesPerProbe * numCorners * sizeof(ProbeSampleData)
    );

    // ------------------------------------------------------------
    // Compute corner metrics normally
    // ------------------------------------------------------------
    std::vector<VoxelLineCornerData> cornerData(numCorners);

    const float3 kLuma = float3(0.2126f, 0.7152f, 0.0722f);

    auto computeEigenvalues3x3Local = [](const float3x3& A, float& e1, float& e2, float& e3)
        {
            const float ONE_THIRD = 1.0f / 3.0f;
            float m = (A[0][0] + A[1][1] + A[2][2]) * ONE_THIRD;
            float k00 = A[0][0] - m, k11 = A[1][1] - m, k22 = A[2][2] - m;
            float k01 = A[0][1], k02 = A[0][2], k12 = A[1][2];
            float q = 0.5f * (k00 * (k11 * k22 - k12 * k12) - k01 * (k01 * k22 - k12 * k02) + k02 * (k01 * k12 - k11 * k02));
            float p = (k00 * k00 + k11 * k11 + k22 * k22 + 2.0f * (k01 * k01 + k02 * k02 + k12 * k12)) / 6.0f;

            if (p < 1e-20f) { e1 = e2 = e3 = m; return; }

            float pSqrt = std::sqrt(p);
            float detVal = q / (p * pSqrt);
            detVal = std::max(-1.0f, std::min(1.0f, detVal));
            float phi = ONE_THIRD * std::acos(detVal);
            float twoSqrtP = 2.0f * pSqrt;
            float s = std::sin(phi), c = std::cos(phi);

            e1 = m + twoSqrtP * c;
            e2 = m - twoSqrtP * (c * 0.5f + s * 0.8660254f);
            e3 = m - twoSqrtP * (c * 0.5f - s * 0.8660254f);
        };

    for (uint32_t cornerIdx = 0; cornerIdx < numCorners; ++cornerIdx)
    {
        int offset = cornerIdx * numSamplesPerProbe;

        std::vector<ProbeSampleData> probeSamplingResults;
        probeSamplingResults.reserve(numSamplesPerProbe);
        for (int s = 0; s < numSamplesPerProbe; ++s)
            probeSamplingResults.push_back(allSamplingData[offset + s]);

        std::vector<float3> coeffs;
        std::vector<GradSHCoeff> grads;
        std::vector<float3x3> lumHessians;

        float3 xPolar;
        xPolar.x = uniqueCornerPositions[cornerIdx].z;
        xPolar.y = uniqueCornerPositions[cornerIdx].x;
        xPolar.z = uniqueCornerPositions[cornerIdx].y;

        calculateSHCoeffs(coeffs, probeSamplingResults, numSamplesPerProbe);
        calculateSHCoeffsGradientsRGBAndHessiansLum(grads, lumHessians, xPolar, probeSamplingResults, samplingDirs);

        VoxelLineCornerData cd{};
        cd.position = uniqueCornerPositions[cornerIdx];

        float sumSqL = 0.0f;
        for (const auto& coeff : coeffs)
        {
            float lum = dot(coeff, kLuma);
            sumSqL += lum * lum;
        }
        cd.coeffVecL2 = std::sqrt(sumSqL);
        cd.isValid = (cd.coeffVecL2 >= 0.0001f);

        float sumSquares = 0.0f;
        for (const auto& h : lumHessians)
        {
            float e1, e2, e3;
            computeEigenvalues3x3Local(h, e1, e2, e3);
            float maxLambda = std::max({ std::abs(e1), std::abs(e2), std::abs(e3) });
            sumSquares += maxLambda * maxLambda;
        }
        cd.maxLambdaVecL2 = std::sqrt(sumSquares);

        cornerData[cornerIdx] = cd;
    }

    // ------------------------------------------------------------
    // Compute voxel error and export CSV
    // ------------------------------------------------------------
    std::ofstream csv("distance_voxel_line_experiment.csv");
    csv << "idx,cx,cy,cz,r,coeffL2_avg,maxLambda_avg,Eabs,Erel\n";

    const float distSq = dot(voxelSize, voxelSize);

    for (uint32_t i = 0; i < (uint32_t)voxels.size(); ++i)
    {
        const auto& voxel = voxels[i];

        float totalCoeff = 0.0f;
        float totalLambda = 0.0f;
        float totalErrorAbs = 0.0f;
        float totalErrorRel = 0.0f;
        uint32_t validCount = 0;

        for (int k = 0; k < 8; ++k)
        {
            const auto& c = cornerData[voxel.cornerIndices[k]];
            if (!c.isValid) continue;

            validCount++;
            totalCoeff += c.coeffVecL2;
            totalLambda += c.maxLambdaVecL2;

            float e = 0.5f * c.maxLambdaVecL2 * distSq;
            totalErrorAbs += e;
            totalErrorRel += e / std::max(c.coeffVecL2, 1e-5f);
        }

        float coeffAvg = validCount > 0 ? totalCoeff / float(validCount) : 0.0f;
        float lambdaAvg = validCount > 0 ? totalLambda / float(validCount) : 0.0f;
        float eAbs = validCount > 0 ? totalErrorAbs / float(validCount) : 0.0f;
        float eRel = validCount > 0 ? totalErrorRel / float(validCount) : 0.0f;

        float r = length(voxel.center - lightPos);

        csv << i << ","
            << voxel.center.x << ","
            << voxel.center.y << ","
            << voxel.center.z << ","
            << r << ","
            << coeffAvg << ","
            << lambdaAvg << ","
            << eAbs << ","
            << eRel << "\n";
    }

    csv.close();
    logInfo("Exported distance_voxel_line_experiment.csv");
}
