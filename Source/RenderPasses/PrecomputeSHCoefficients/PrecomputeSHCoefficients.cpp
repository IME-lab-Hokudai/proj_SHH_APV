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
#include "Scene/TriangleMesh.h"

const int sampleCount = 4096;

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
    reflector.addOutput("output", "Color").texture2D(sz.x, sz.y, 4);
    // Add the required depth output. This always exists.
    reflector.addOutput("depth", "Depth buffer")
        .format(ResourceFormat::D32Float)
        .bindFlags(ResourceBindFlags::DepthStencil)
        .texture2D(sz.x, sz.y, 4);
    //reflector.addInput()
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

        pRenderContext->clearDsv(pDepth->getDSV().get(), 1.f, 0);

        auto shShaderRootVar = mpVars->getRootVar();
        shShaderRootVar["gLinearSampler"] = mpLinearSampler;
        //shShaderRootVar["PerFrameCB"]["shCoeffs"].setBlob(shCoeffs.data(), shCoeffs.size() * sizeof(float4)); // bind sh coeffs to cbuffer
        shShaderRootVar["gSHCoeffs"] = mpGridSHBuffer;
        shShaderRootVar["gProbeGridInfo"]["resolution"] = mProbeGrid.resolution;
        shShaderRootVar["gProbeGridInfo"]["numBasis"] = mProbeGrid.numBasis;
        shShaderRootVar["gProbeGridInfo"]["origin"] = mProbeGrid.origin;
        shShaderRootVar["gProbeGridInfo"]["spacing"] = mProbeGrid.spacing;

        mpScene->rasterize(pRenderContext, mpGraphicsState.get(), mpVars.get(), mpRasterState, mpRasterState);

        if (!mbFinishProbeSampling)
        {
            auto rtVar = mpRtVars->getRootVar();
            rtVar["gProbeDirSamples"] = mpProbeDirSamplesBuffer;
            rtVar["gProbeSamplingOutput"] = mpProbeSamplingResultBuffer;
            rtVar["PerFrameCB"]["sampleCount"] = sampleCount;
            rtVar["PerFrameCB"]["sampleIndex"] = mSampleIndex++;
            mpScene->raytrace(pRenderContext, mpRtProgram.get(), mpRtVars, uint3(sampleCount, 1, 1));

            // Ensure all UAV writes are finished before mapping for read-back
           // pRenderContext->uavBarrier(mpProbeSamplingResultBuffer.get());
           // mpDevice->wait();

            // Map the buffer for reading
            //float4* pData = new float4[sampleCount];
            //mpProbeSamplingResultBuffer->getBlob(pData, 0, sampleCount*sizeof(float4));

            //// Now pData points to your results, size is sampleCount
            //for (int i = 0; i < sampleCount; ++i)
            //{
            //    float4 result = pData[i];
            //    // Process result as needed
            //    logInfo(fmt::format("AAAA: {:.3f} {:.3f} {:.3f}", result.x, result.y, result.z));
            //}

            // Unmap when done
            mbFinishProbeSampling = true;
           // delete[] pData;
        }

        if (mbFinishProbeSampling)
        {
            // if (mbShowSHGrid)
            //{
            mpProbeVisualizePass->setCameraData(
                mpScene->getCamera()->getViewProjMatrix(), mpScene->getCamera()->getViewMatrix(), mpScene->getCamera()->getProjMatrix()
            );
            mpProbeVisualizePass->setProbeSamplingData(mpProbeDirSamplesBuffer, mpProbeSamplingResultBuffer);
            mpProbeVisualizePass->execute(pRenderContext, mpFbo);
            //}
        }

    }
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
      //if (mpScene->getEnvMap() != nullptr)
        //{
            mpEnvMap = mpScene->getEnvMap();
            initSHTable(2, mpEnvMap->getEnvMap()->getWidth(), mpEnvMap->getEnvMap()->getHeight());
            decomposeSH(shCoeffs, mpEnvMap);
            //reconstructSH(shCoeffs, mpEnvMap, mpDevice);

            //mProbeGrid.origin = float3(0.0f, 0.0f, 0.0f);
            //mProbeGrid.resolution = int3(3, 3, 3);
            //mProbeGrid.spacing = float3(0.5f, 0.5f, 0.5f);

            //createProbeGrid(mProbeGrid, shCoeffs);
            //saveProbeGridToFile(mProbeGrid, "ProbeGrid.txt");

            //loadProbeGridFromFile(mProbeGrid, "ProbeGrid.txt");

            AABB sceneBounds = mpScene->getSceneBounds();
            float3 minBound = sceneBounds.minPoint;
            float3 maxBound = sceneBounds.maxPoint;
            float3 sceneCenter = sceneBounds.center();
            float3 sceneSize = maxBound - minBound;

            // Decide spacing between probes
            float3 spacing = float3(1.0f, 1.0f, 1.0f); 
            // Number of probes in each dimension
            int3 resolution;
            resolution.x = (int)ceil(sceneSize.x / spacing.x) + 1;
            resolution.y = (int)ceil(sceneSize.y / spacing.y) + 1;
            resolution.z = (int)ceil(sceneSize.z / spacing.z) + 1;

            float3 halfSize = 0.5f*(float3(resolution) - 1.0f) * spacing;
            resolution = int3(1, 1, 1);
            //mProbeGrid.origin = sceneCenter;
            //mProbeGrid.origin = sceneCenter - halfSize;
            mProbeGrid.origin = sceneCenter;
            mProbeGrid.spacing = spacing;
            mProbeGrid.resolution = resolution;

            createProbeGrid(mProbeGrid, shCoeffs);

           // shCoeffs.clear();
           // shCoeffs.insert(shCoeffs.end(), mProbeGrid.probes.begin(), mProbeGrid.probes.begin() + mProbeGrid.numBasis); //just one env set for envmap render

           int numProbes = mProbeGrid.resolution.x * mProbeGrid.resolution.y * mProbeGrid.resolution.z;
            mpGridSHBuffer = mpDevice->createStructuredBuffer(
                sizeof(float4),
                mProbeGrid.numBasis * numProbes,
                ResourceBindFlags::ShaderResource,
                MemoryType::DeviceLocal,
               mProbeGrid.probes.data()
            );
           mpGridSHBuffer->setName("SH Grid Info");

            // float3 probePos = sceneCenter - halfSize;
            auto dirSamples = generateUniformSphereDirSamples(sampleCount, sceneCenter);

            mpProbeDirSamplesBuffer = mpDevice->createStructuredBuffer(
                sizeof(ProbeDirSample), sampleCount, ResourceBindFlags::ShaderResource, MemoryType::DeviceLocal, dirSamples.data()
            );
            mpProbeDirSamplesBuffer->setName("Probe Dir Samples");
            mpProbeSamplingResultBuffer = mpDevice->createStructuredBuffer(
                sizeof(float4), sampleCount, ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess, MemoryType::DeviceLocal
            );
            mpProbeSamplingResultBuffer->setName("Probe Sampling Result Buffer");

           //mpFullScreenPass = FullScreenPass::create(mpDevice, kEnvMapShaderFile, mpScene->getSceneDefines(), 0, "vsMain");
           mpProbeVisualizePass = ProbeVisualizePass::create(mpDevice, mpScene->getSceneDefines());
           mpProbeVisualizePass->setGridData(mProbeGrid, dirSamples);

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

        // ray tracing program to sample probes
        ProgramDesc rtProgDesc;
        rtProgDesc.addShaderModules(mpScene->getShaderModules());
        rtProgDesc.addShaderLibrary(kProbeSamplingFile);
        rtProgDesc.setMaxTraceRecursionDepth(3); // 1 for calling TraceRay from RayGen, 1 for calling it from the
                                                 // primary-ray ClosestHit shader for reflections, 1 for reflection ray
                                                 // tracing a shadow ray
        rtProgDesc.setMaxPayloadSize(24);        // The largest ray payload struct (PrimaryRayData) is 24 bytes. The payload size
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

        mpRtProgram = Program::create(mpDevice, rtProgDesc, mpScene->getSceneDefines());
        mpRtVars = RtProgramVars::create(mpDevice, mpRtProgram, sbt);
    }
}
