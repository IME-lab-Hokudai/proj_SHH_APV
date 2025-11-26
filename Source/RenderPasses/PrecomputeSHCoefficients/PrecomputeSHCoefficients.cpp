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
#include "PrecomputeSHCoefficients.h"

#include "envMap_SH.h"
#include "Core/Pass/FullScreenPass.h"
#include "Rendering/Lights/EmissivePowerSampler.h"
#include "Rendering/Lights/EmissiveUniformSampler.h"
#include "Rendering/Lights/LightBVHSampler.h"
#include "Scene/SceneBuilder.h"
#include "Scene/TriangleMesh.h"
#include "ProbeSamplingData.slang"
const int numSamplesPerProbe = 4096;
const int verificationRes = 32;
const float verificationH = 0.005f;
const float verificationY = 0.2f;
const float verificationExtent = 0.25f;
namespace
{
//const char kShaderFile[] = "RenderPasses/PrecomputeSHCoefficients/SHShader.slang";
const char kShaderFile[] = "RenderPasses/PrecomputeSHCoefficients/SHGridShader.slang";
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
            mbShowSHGrid = value;
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
    props[kShowSHGrid] = mbShowSHGrid;
    return props;
}

RenderPassReflection PrecomputeSHCoefficients::reflect(const CompileData& compileData)
{
    // Define the required resources here
    RenderPassReflection reflector;
    const uint2 sz = RenderPassHelpers::calculateIOSize(mOutputSizeSelection, mFixedOutputSize, compileData.defaultTexDims);
    // REMARK MSAA is set via texture sample count. Note that all fbo attachment have to have same sample count.
    //reflector.addOutput("output", "Color").texture2D(sz.x, sz.y, 4).format(ResourceFormat::RGBA16Float);
    // Add the required depth output. This always exists.
    // reflector.addOutput("depth", "Depth buffer")
    //     .format(ResourceFormat::D32Float)
    //     .bindFlags(ResourceBindFlags::DepthStencil)
    //     .texture2D(sz.x, sz.y, 4);

    reflector.addOutput("output", "Color").texture2D(sz.x, sz.y).format(ResourceFormat::RGBA16Float);
    reflector.addOutput("depth", "Depth buffer")
        .format(ResourceFormat::D32Float)
        .bindFlags(ResourceBindFlags::DepthStencil)
        .texture2D(sz.x, sz.y);
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
        //auto envMapShaderRootVar = mpFullScreenPass->getRootVar();

        //mpEnvMap->bindShaderData(envMapShaderRootVar["gScene"]["envMap"]);
        //mpScene->getCamera()->bindShaderData(envMapShaderRootVar["gScene"]["camera"]);
        //envMapShaderRootVar["PerFrameCB"]["shCoeffs"].setBlob(shCoeffs.data(), shCoeffs.size() * sizeof(float4)); // bind sh coeffs to cbuffer
        //envMapShaderRootVar["PerFrameCB"]["showReconstructedEnvMap"] = mbShowReconstructedEnvMap;
        //mpFullScreenPass->execute(pRenderContext, mpFbo);

 //       if (mbVerify)
 //       {
 //           auto rtVar = mpRtVars->getRootVar();
 //           rtVar["gProbeDirSamples"] = mpProbeDirSamplesBuffer;
 //           rtVar["gProbePositions"] = mpProbePosBuffer;

 //           if (mpEmissiveSampler)
 //               mpEmissiveSampler->bindShaderData(rtVar["PerFrameCB"]["emissiveSampler"]);

 //           rtVar["gProbeSamplingOutput"] = mpProbeSamplingResultBuffer;
 //           rtVar["PerFrameCB"]["numSamplePerProbe"] = numSamplesPerProbe;
 //           rtVar["PerFrameCB"]["sampleIndex"] = mSampleIndex++;
 //           //int numProbe = verificationRes*verificationRes*3;
 //           int numProbe = verificationRes * 3; // one scanline only
 //           mpScene->raytrace(pRenderContext, mpRtProgram.get(), mpRtVars, uint3(numSamplesPerProbe, numProbe, 1));

 //           ProbeSampleData* allProbeSamplingData = new ProbeSampleData[numSamplesPerProbe * numProbe];
 //           mpProbeSamplingResultBuffer->getBlob(allProbeSamplingData, 0, numSamplesPerProbe * numProbe * sizeof(ProbeSampleData));

 //           int basisIdx = 0; // l=2, m=0 l(l+1)+m = 2*3 + 0 = 6
 //           int numBasis = 9;
 //           std::vector<float> coeffVec;
 //           coeffVec.clear();
 //           coeffVec.reserve(numProbe);

 //            std::vector<float> finiteDifferenceGrads;
 //            finiteDifferenceGrads.clear();
 //            finiteDifferenceGrads.reserve(numProbe/3);

 //            std::vector<float> finiteDifferenceHessians;
 //            finiteDifferenceHessians.clear();
 //            finiteDifferenceHessians.reserve(numProbe / 3);

 //             std::vector<float3> analyticGrads;
 //             analyticGrads.clear();
 //             analyticGrads.reserve(numProbe / 3);

 //             std::vector<float> analyticGradsY;
 //             analyticGradsY.clear();
 //             analyticGradsY.reserve(numProbe / 3);

 //             std::vector<float> analyticHessYY;
 //             analyticHessYY.clear();
 //             analyticHessYY.reserve(numProbe / 3);

 //           float verificationHSq = verificationH * verificationH;
 //           // calculate SH coeffs for all verification positions
 //           for (int probeIdx = 0; probeIdx < numProbe; probeIdx+=3)
 //           {
 //               int offset = probeIdx * numSamplesPerProbe;
 //               std::vector<ProbeSampleData> probeSamplingResults;
 //               probeSamplingResults.clear();
 //               probeSamplingResults.reserve(numSamplesPerProbe);
 //               for (int sampleIdx = 0; sampleIdx < numSamplesPerProbe; sampleIdx++)
 //               {
 //                   probeSamplingResults.push_back(allProbeSamplingData[offset + sampleIdx]);
 //               }
 //               float coeffR = calculateChannelRSHCoeffLM(basisIdx, probeSamplingResults);
 //               coeffVec.push_back(coeffR);

 //               float3 xPolar; // x in polar coord
 //               xPolar.x = verificationPositions[probeIdx].z;
 //               xPolar.y = verificationPositions[probeIdx].x;
 //               xPolar.z = verificationPositions[probeIdx].y;

 //               float3 xPolarPPlus; // x' in polar coord
 //               xPolarPPlus.x = verificationPositions[probeIdx + 1].z;
 //               xPolarPPlus.y = verificationPositions[probeIdx + 1].x;
 //               xPolarPPlus.z = verificationPositions[probeIdx + 1].y;

 //               float3 xPolarPMinus; // x' in polar coord
 //               xPolarPMinus.x = verificationPositions[probeIdx + 2].z;
 //               xPolarPMinus.y = verificationPositions[probeIdx + 2].x;
 //               xPolarPMinus.z = verificationPositions[probeIdx + 2].y;
 //               
 //               float coeffRPPlus = calculateCoeffRPrime(probeSamplingResults, xPolar, xPolarPPlus, numBasis, basisIdx);
 //               float coeffRPMinus = calculateCoeffRPrime(probeSamplingResults, xPolar, xPolarPMinus, numBasis, basisIdx);

 //               float finiteGrad = (coeffRPPlus - coeffR) / verificationH; // ∂f/∂x≈(f(x+h)-f(x))/h
 //               //float finiteGrad = (coeffR - coeffRPMinus) / verificationH;
 //               finiteDifferenceGrads.push_back(finiteGrad);

 //               float finiteHessYY = (coeffRPPlus - 2.0f * coeffR + coeffRPMinus)/verificationHSq;
 //               finiteDifferenceHessians.push_back(finiteHessYY);

 //               float3 analyticGrad = float3(0.0f, 0.0f, 0.0f);
 //               float3x3 analyticHess = float3x3::zeros();
 //               calculateChannelRGradAndHessianSHCoeffLM(xPolar, probeSamplingResults, samplingDirs, basisIdx, analyticGrad, analyticHess);
 //               //computeKrivanekCoeffLMGradient(xPolar, probeSamplingResults, basisIdx, analyticGrad);
 //               analyticGradsY.push_back(analyticGrad.y);
 //               analyticHessYY.push_back(analyticHess[1][1]);
 //           }

 ///*           for (int i = 0; i < analyticGrads.size(); i++)
 //           {
 //                analyticGradY.push_back(analyticGrads[i].y);
 //           }*/


 //           mbVerify = false;
 //           delete[] allProbeSamplingData;
 //       }

        if (mbVerify) // mixed derivative verification
        {
            auto rtVar = mpRtVars->getRootVar();
            rtVar["gProbeDirSamples"] = mpProbeDirSamplesBuffer;
            rtVar["gProbePositions"] = mpProbePosBuffer;

            if (mpEmissiveSampler)
                mpEmissiveSampler->bindShaderData(rtVar["PerFrameCB"]["emissiveSampler"]);

            rtVar["gProbeSamplingOutput"] = mpProbeSamplingResultBuffer;
            rtVar["PerFrameCB"]["numSamplePerProbe"] = numSamplesPerProbe;
            rtVar["PerFrameCB"]["sampleIndex"] = mSampleIndex++;
            // int numProbe = verificationRes*verificationRes*3
            int numProbe = verificationRes * 5; // one scanline only
            mpScene->raytrace(pRenderContext, mpRtProgram.get(), mpRtVars, uint3(numSamplesPerProbe, numProbe, 1));

            ProbeSampleData* allProbeSamplingData = new ProbeSampleData[numSamplesPerProbe * numProbe];
            mpProbeSamplingResultBuffer->getBlob(allProbeSamplingData, 0, numSamplesPerProbe * numProbe * sizeof(ProbeSampleData));

            int basisIdx = 2; // l=2, m=0 l(l+1)+m = 2*3 + 0 = 6
            int numBasis = 9;
            std::vector<float> coeffVec;
            coeffVec.clear();
            coeffVec.reserve(numProbe);

            std::vector<float> finiteDifferenceGrads;
            finiteDifferenceGrads.clear();
            finiteDifferenceGrads.reserve(numProbe / 3);

            std::vector<float> finiteDifferenceHessians;
            finiteDifferenceHessians.clear();
            finiteDifferenceHessians.reserve(numProbe / 3);

            std::vector<float3> analyticGrads;
            analyticGrads.clear();
            analyticGrads.reserve(numProbe / 3);

            std::vector<float> analyticHessYX;
            analyticHessYX.clear();
            analyticHessYX.reserve(numProbe / 3);

            float verificationHSq =4.0f * verificationH * verificationH;
            // calculate SH coeffs for all verification positions
            for (int probeIdx = 0; probeIdx < numProbe; probeIdx += 5)
            {
                int offset = probeIdx * numSamplesPerProbe;
                std::vector<ProbeSampleData> probeSamplingResults;
                probeSamplingResults.clear();
                probeSamplingResults.reserve(numSamplesPerProbe);
                for (int sampleIdx = 0; sampleIdx < numSamplesPerProbe; sampleIdx++)
                {
                    probeSamplingResults.push_back(allProbeSamplingData[offset + sampleIdx]);
                }
                //float coeffR = calculateChannelRSHCoeffLM(basisIdx, probeSamplingResults);
                //coeffVec.push_back(coeffR);

                float3 xPolar; // x in polar coord
                xPolar.x = verificationPositions[probeIdx].z;
                xPolar.y = verificationPositions[probeIdx].x;
                xPolar.z = verificationPositions[probeIdx].y;

                float3 xPolarXPlusZPlus; // x' in polar coord
                xPolarXPlusZPlus.x = verificationPositions[probeIdx + 1].z;
                xPolarXPlusZPlus.y = verificationPositions[probeIdx + 1].x;
                xPolarXPlusZPlus.z = verificationPositions[probeIdx + 1].y;

                float3 xPolarXPlusZMinus; // x' in polar coord
                xPolarXPlusZMinus.x = verificationPositions[probeIdx + 2].z;
                xPolarXPlusZMinus.y = verificationPositions[probeIdx + 2].x;
                xPolarXPlusZMinus.z = verificationPositions[probeIdx + 2].y;

                float3 xPolarXMinusZPlus; // x' in polar coord
                xPolarXMinusZPlus.x = verificationPositions[probeIdx + 3].z;
                xPolarXMinusZPlus.y = verificationPositions[probeIdx + 3].x;
                xPolarXMinusZPlus.z = verificationPositions[probeIdx + 3].y;

                 float3 xPolarXMinusZMinus; // x' in polar coord
                xPolarXMinusZMinus.x = verificationPositions[probeIdx + 4].z;
                 xPolarXMinusZMinus.y = verificationPositions[probeIdx + 4].x;
                xPolarXMinusZMinus.z = verificationPositions[probeIdx + 4].y;

                float coeffRXPlusZPlus = calculateCoeffRPrime(probeSamplingResults, xPolar, xPolarXPlusZPlus, numBasis, basisIdx);
                float coeffRXPlusZMinus = calculateCoeffRPrime(probeSamplingResults, xPolar, xPolarXPlusZMinus, numBasis, basisIdx);
                float coeffRXMinusZPlus = calculateCoeffRPrime(probeSamplingResults, xPolar, xPolarXMinusZPlus, numBasis, basisIdx);
                float coeffRXMinusZMinus = calculateCoeffRPrime(probeSamplingResults, xPolar, xPolarXMinusZMinus, numBasis, basisIdx);

                float finiteHessYX = (coeffRXPlusZPlus - coeffRXPlusZMinus - coeffRXMinusZPlus + coeffRXMinusZMinus) / verificationHSq;
                finiteDifferenceHessians.push_back(finiteHessYX);

                float3 analyticGrad = float3(0.0f, 0.0f, 0.0f);
                float3x3 analyticHess = float3x3::zeros();
                calculateChannelRGradAndHessianSHCoeffLM(xPolar, probeSamplingResults, samplingDirs, basisIdx, analyticGrad, analyticHess);
                analyticHessYX.push_back(analyticHess[0][1]);
            }

            mbVerify = false;
            delete[] allProbeSamplingData;
        }

        if (!mbFinishSHPrecompute)
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

        // visualize probes
        if (mbFinishSHPrecompute)
        {
            //pRenderContext->clearDsv(pDepth->getDSV().get(), 1.f, 0);

            // auto shShaderRootVar = mpVars->getRootVar();
            // shShaderRootVar["gLinearSampler"] = mpLinearSampler;
            // shShaderRootVar["gSHCoeffs"] = mpGridSHCoeffsBuffer;
            // shShaderRootVar["gProbeGridInfo"]["resolution"] = mProbeGrid.resolution;
            // shShaderRootVar["gProbeGridInfo"]["numBasis"] = mProbeGrid.numBasis;
            // shShaderRootVar["gProbeGridInfo"]["origin"] = mProbeGrid.origin;
            // shShaderRootVar["gProbeGridInfo"]["spacing"] = mProbeGrid.spacing;

            // mpScene->rasterize(pRenderContext, mpGraphicsState.get(), mpVars.get(), mpRasterState, mpRasterState);

            // TODO: rework grid visualization pass
             //if (mbShowSHGrid)
             //{
             //   mpProbeVisualizePass->setCameraData(
             //       mpScene->getCamera()->getViewProjMatrix(), mpScene->getCamera()->getViewMatrix(), mpScene->getCamera()->getProjMatrix()
             //   );
             //   mpProbeVisualizePass->setProbeSamplingData(mpProbeDirSamplesBuffer, mpProbeSamplingResultBuffer);
             //   mpProbeVisualizePass->execute(pRenderContext, mpFbo);
             //}
        }
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
    float Omega_i = (4.0f * M_PI) / (float)numSamplesPerProbe;

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
    if (widget.checkbox("Show SH Grid", mbShowSHGrid))
        requestRecompile();
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

        //generate verification data
           if (mbVerify)
           {
               //std::vector<float3> verificationPositions;
               if (sceneName == "cornell")
               {
                  // AABB sceneBounds = mpScene->getSceneBounds();
                  // float3 minBound = sceneBounds.minPoint;
                  // float3 maxBound = sceneBounds.maxPoint;
                  // float3 sceneCenter = sceneBounds.center();
                  // float3 sceneSize = maxBound - minBound;

                  //float3 spacing = float3(.15f, .15f, .15f);
                  //// Number of probes in each dimension
                  // int3 resolution;
                  // resolution.x = (int)ceil(sceneSize.x / spacing.x);
                  // resolution.y = (int)ceil(sceneSize.y / spacing.y);
                  // resolution.z = (int)ceil(sceneSize.z / spacing.z);

                  //// resolution = int3(1, 1, 1); // for testing

                  // float3 halfSize = 0.5f * (float3(resolution) - 1.0f) * spacing;
                  // // resolution = int3(1, 1, 1);
                  // mProbeGrid.origin = sceneCenter - halfSize;
                  // mProbeGrid.origin += 0.025f;

                  // // mProbeGrid.origin = sceneCenter;
                  // mProbeGrid.spacing = spacing;
                  // mProbeGrid.resolution = resolution;
                  // mProbeGrid.numBasis = 9; // number of SH coefficients per probe
                  // computeProbesPos(mProbeGrid);

                  // 

                  // verificationPositions.clear();
                  // verificationPositions.reserve(3);

                  // float3 tmpX = mProbeGrid.probesPos[0];
                  // float3 tmpXForward = tmpX + float3(verificationH, 0.0f, 0.0f);
                  // float3 tmpXBackward = tmpX - float3(verificationH, 0.0f, 0.0f);
                  // verificationPositions.push_back(tmpX);
                  // verificationPositions.push_back(tmpXForward);
                  // verificationPositions.push_back(tmpXBackward);

                   //verificationPositions = generateVerificationPositions(verificationY, verificationExtent, verificationRes, verificationH); // Falcor coord sample: [center,+h, -h]
                   verificationPositions = generateVerificationPositionsMixed(verificationY, verificationExtent, verificationRes, verificationH); // falcor coord Order per sample: [center, x+h z+h, x+h z-h, x-h z+h, x-h z-h]
                   int numProbes = verificationPositions.size();
                   mpProbePosBuffer = mpDevice->createStructuredBuffer(
                       sizeof(float3), numProbes, ResourceBindFlags::ShaderResource, MemoryType::DeviceLocal, verificationPositions.data()
                   );
                   mpProbePosBuffer->setName("probes world pos");

                    mpProbeSamplingResultBuffer = mpDevice->createStructuredBuffer(
                    sizeof(ProbeSampleData),
                       numSamplesPerProbe * numProbes,
                    ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess,
                    MemoryType::DeviceLocal
                    );
                    mpProbeSamplingResultBuffer->setName("Probe Sampling Result Buffer");
               }
               ProgramDesc rtProgDesc;
               rtProgDesc.addShaderModules(mpScene->getShaderModules());
               rtProgDesc.addShaderLibrary(kProbeSamplingFile);
               rtProgDesc.setMaxTraceRecursionDepth(3); // 1 for calling TraceRay from RayGen, 1 for calling it from the
                                                        // primary-ray ClosestHit shader for reflections, 1 for reflection ray
                                                        // tracing a shadow ray
               rtProgDesc.setMaxPayloadSize(64);        // The largest ray payload struct (PrimaryRayData) is 24 bytes. The payload size
                                                        // should be set as small as possible for maximum performance.
               rtProgDesc.setMaxAttributeSize(8);
               // Add global type conformances.
               rtProgDesc.addTypeConformances(mpScene->getTypeConformances());

               ref<RtBindingTable> sbt = RtBindingTable::create(2, 2, mpScene->getGeometryCount());
               sbt->setRayGen(rtProgDesc.addRayGen("rayGen"));
               sbt->setMiss(0, rtProgDesc.addMiss("primaryMiss"));
               // sbt->setMiss(1, rtProgDesc.addMiss("shadowMiss"));
               auto primary = rtProgDesc.addHitGroup("primaryClosestHit");
               // auto shadow = rtProgDesc.addHitGroup("", "shadowAnyHit");

               sbt->setHitGroup(0, mpScene->getGeometryIDs(Scene::GeometryType::TriangleMesh), primary);
               //  sbt->setHitGroup(1, mpScene->getGeometryIDs(Scene::GeometryType::TriangleMesh), shadow);

               const auto& pLights = mpScene->getILightCollection(pRenderContext); // REMARK weird design that light collection is created
                                                                                   // upon first call to this.
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
           }

           //if (!mbFinishSHPrecompute) // config for precomputing SH coefficients
           //{
           //    AABB sceneBounds = mpScene->getSceneBounds();
           //    float3 minBound = sceneBounds.minPoint;
           //    float3 maxBound = sceneBounds.maxPoint;
           //    float3 sceneCenter = sceneBounds.center();
           //    float3 sceneSize = maxBound - minBound;
           //   
           //    int order = 2; // SH order

           //    initSHBasisGradientAndHessianTables(dirSamples);

           //    float3 spacing = float3(1.f, 1.f, 1.f);

           //    if (sceneName == "arcade")
           //    {
           //        spacing = float3(1.f, 1.f, 1.f);
           //     }
           //     else
           //     {
           //         spacing = float3(.15f, .15f, .15f);
           //    }
           //    // Number of probes in each dimension
           //    int3 resolution;
           //    resolution.x = (int)ceil(sceneSize.x / spacing.x);
           //    resolution.y = (int)ceil(sceneSize.y / spacing.y);
           //    resolution.z = (int)ceil(sceneSize.z / spacing.z);

           //   // resolution = int3(1, 1, 1); // for testing

           //    float3 halfSize = 0.5f * (float3(resolution) - 1.0f) * spacing;
           //    // resolution = int3(1, 1, 1);
           //    mProbeGrid.origin = sceneCenter - halfSize;
           //    mProbeGrid.origin += 0.025f;

           //    // mProbeGrid.origin = sceneCenter;
           //    mProbeGrid.spacing = spacing;
           //    mProbeGrid.resolution = resolution;
           //    mProbeGrid.numBasis = (order + 1) * (order + 1); // number of SH coefficients per probe
           //    int numProbes = mProbeGrid.resolution.x * mProbeGrid.resolution.y * mProbeGrid.resolution.z;
           //    computeProbesPos(mProbeGrid);

           //    mpProbeSamplingResultBuffer = mpDevice->createStructuredBuffer(
           //        sizeof(ProbeSampleData),
           //        numSamplePerProbe * numProbes,
           //        ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess,
           //        MemoryType::DeviceLocal
           //    );
           //    mpProbeSamplingResultBuffer->setName("Probe Sampling Result Buffer");

           //    mpProbePosBuffer = mpDevice->createStructuredBuffer(
           //        sizeof(float3), numProbes, ResourceBindFlags::ShaderResource, MemoryType::DeviceLocal, mProbeGrid.probesPos.data()
           //    );
           //    mpProbePosBuffer->setName("probes world pos");

           //            // ray tracing program to sample probes
           //    ProgramDesc rtProgDesc;
           //    rtProgDesc.addShaderModules(mpScene->getShaderModules());
           //    rtProgDesc.addShaderLibrary(kProbeSamplingFile);
           //    rtProgDesc.setMaxTraceRecursionDepth(3); // 1 for calling TraceRay from RayGen, 1 for calling it from the
           //                                             // primary-ray ClosestHit shader for reflections, 1 for reflection ray
           //                                             // tracing a shadow ray
           //    rtProgDesc.setMaxPayloadSize(64);        // The largest ray payload struct (PrimaryRayData) is 24 bytes. The payload size
           //                                             // should be set as small as possible for maximum performance.
           //    rtProgDesc.setMaxAttributeSize(8);
           //    // Add global type conformances.
           //    rtProgDesc.addTypeConformances(mpScene->getTypeConformances());

           //    ref<RtBindingTable> sbt = RtBindingTable::create(2, 2, mpScene->getGeometryCount());
           //    sbt->setRayGen(rtProgDesc.addRayGen("rayGen"));
           //    sbt->setMiss(0, rtProgDesc.addMiss("primaryMiss"));
           //    // sbt->setMiss(1, rtProgDesc.addMiss("shadowMiss"));
           //    auto primary = rtProgDesc.addHitGroup("primaryClosestHit");
           //    // auto shadow = rtProgDesc.addHitGroup("", "shadowAnyHit");

           //    sbt->setHitGroup(0, mpScene->getGeometryIDs(Scene::GeometryType::TriangleMesh), primary);
           //    //  sbt->setHitGroup(1, mpScene->getGeometryIDs(Scene::GeometryType::TriangleMesh), shadow);

           //    const auto& pLights = mpScene->getILightCollection(pRenderContext); //REMARK weird design that light collection is created upon first call to this.
           //    if (mpScene->useEmissiveLights())
           //    {
           //        if (!mpEmissiveSampler)
           //        {
           //            FALCOR_ASSERT(pLights && pLights->getActiveLightCount(pRenderContext) > 0);
           //            FALCOR_ASSERT(!mpEmissiveSampler);

           //            switch (mEmissiveSamplerType)
           //            {
           //                case EmissiveLightSamplerType::Uniform: // use uniform sampling as default for now
           //                    mpEmissiveSampler =
           //                        std::make_unique<EmissiveUniformSampler>(pRenderContext, mpScene->getILightCollection(pRenderContext));
           //                    break;
           //                case EmissiveLightSamplerType::LightBVH:
           //                    mpEmissiveSampler = std::make_unique<LightBVHSampler>(
           //                        pRenderContext, mpScene->getILightCollection(pRenderContext), mLightBVHOptions
           //                    );
           //                    break;
           //                case EmissiveLightSamplerType::Power:
           //                    mpEmissiveSampler =
           //                        std::make_unique<EmissivePowerSampler>(pRenderContext, mpScene->getILightCollection(pRenderContext));
           //                    break;
           //                default:
           //                    FALCOR_THROW("Unknown emissive light sampler type");
           //            }
           //        }
           //    }

           //   mpRtProgram = Program::create(mpDevice, rtProgDesc, mpScene->getSceneDefines());
           //   
           //    if (mpEmissiveSampler)
           //    {
           //        auto defines = mpEmissiveSampler->getDefines();
           //        mpRtProgram->addDefines(defines);
           //    }

           //    DefineList lightRelatedDefines;
           //    lightRelatedDefines.add("USE_ANALYTIC_LIGHTS", mpScene->useAnalyticLights() ? "1" : "0");
           //    lightRelatedDefines.add("USE_EMISSIVE_LIGHTS", mpScene->useEmissiveLights() ? "1" : "0");

           //    mpRtProgram->addDefines(lightRelatedDefines);

           //    mpRtVars = RtProgramVars::create(mpDevice, mpRtProgram, sbt);
           //}
           //else // config to render scene with resulting SH grid stored in file
           //{
           //    //loadProbeGridFromFile(mProbeGrid, "ProbeGrid.txt");
           //    //loadProbeGridFromFile(mProbeGrid, "ProbeGridCornell.txt");
           //    //if (sceneName == "arcade")
           //    //{
           //    //    loadProbeGridFromFile(mProbeGrid, "ProbeGridArcade.txt");
           //    //}
           //    //else
           //    //{
           //    //    loadProbeGridFromFile(mProbeGrid, "ProbeGridCornell.txt");
           //    //}

           //    if (sceneName == "arcade")
           //    {
           //        loadProbeGridFromFileWithGradAndHessian(mProbeGrid, "ProbeGridWithGradAndHessianArcade.txt");
           //    }
           //    else
           //    {
           //        loadProbeGridFromFileWithGradAndHessian(mProbeGrid, "ProbeGridWithGradAndHessianCornell.txt");
           //    }

           //    //int order = (int)sqrt(mProbeGrid.numBasis) - 1; // SH order
           //    //initSHTable(order, dirSamples);
           //    //initSHBasisGradientAndHessianTables(dirSamples);
           //    int numProbes = mProbeGrid.resolution.x * mProbeGrid.resolution.y * mProbeGrid.resolution.z;
           //    mpGridSHCoeffsBuffer = mpDevice->createStructuredBuffer(
           //        sizeof(float4),
           //        mProbeGrid.numBasis * numProbes,
           //        ResourceBindFlags::ShaderResource,
           //        MemoryType::DeviceLocal,
           //        mProbeGrid.probesSHCoeffs.data()
           //    );
           //    mpGridSHCoeffsBuffer->setName("SH Grid Coeffs");

           //    //std::vector<float4> reconstructedData;
           //    //reconstructSH(mProbeGrid, numSamplePerProbe, reconstructedData);

           //    //mpProbeSamplingResultBuffer = mpDevice->createStructuredBuffer(
           //    //    sizeof(float4),
           //    //    numSamplePerProbe * numProbes,
           //    //    ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess,
           //    //    MemoryType::DeviceLocal,
           //    //    reconstructedData.data()
           //    //);
           //    //mpProbeSamplingResultBuffer->setName("Probe Sampling Result Buffer");
           //}
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
           rasterDesc.setDepthBias(100000, 1.0f);
           mpRasterState = RasterizerState::create(rasterDesc);

           // default depth stencil state
           DepthStencilState::Desc dsDesc;
           ref<DepthStencilState> pDsState = DepthStencilState::create(dsDesc);

           mpGraphicsState = GraphicsState::create(mpDevice);
           mpGraphicsState->setProgram(mpProgram);
           mpGraphicsState->setRasterizerState(mpRasterState);
           mpGraphicsState->setFbo(mpFbo);
           mpGraphicsState->setDepthStencilState(pDsState);
            // mpFullScreenPass = FullScreenPass::create(mpDevice, kEnvMapShaderFile, mpScene->getSceneDefines(), 0, "vsMain");
           //mpProbeVisualizePass = ProbeVisualizePass::create(mpDevice, mpScene->getSceneDefines());
           //mpProbeVisualizePass->setGridData(mProbeGrid, dirSamples);
    }
}


