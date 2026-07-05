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
#include <fstream>
#include "AdaptiveSHDemo.h"
#include "Rendering/Lights/EmissivePowerSampler.h"
#include "Rendering/Lights/EmissiveUniformSampler.h"
#include "Utils/Logger.h"

#define PROBE_MODE_ADAPTIVE 0
#define PROBE_MODE_UNIFORM  1
 // CHANGE THIS LINE TO SWITCH MODES:
//#define CURRENT_PROBE_MODE PROBE_MODE_UNIFORM
#define CURRENT_PROBE_MODE PROBE_MODE_ADAPTIVE

#if CURRENT_PROBE_MODE == PROBE_MODE_ADAPTIVE
//const std::string loadFromFileName = "Seeded8DirectAbsErr5SubwayCorridorNoOpen.txt";
//const std::string loadFromFileName = "DirectAbsErr20DataScene.txt";
//const std::string loadFromFileName = "DirectAbsErr10DataScene.txt";
//const std::string loadFromFileName = "DirectAbsErr8p5N6DataScene.txt";
//const std::string loadFromFileName = "DirectAbsErr5DataScene.txt";
//const std::string loadFromFileName = "DirectAbsErr1DataScene.txt";
//const std::string loadFromFileName = "DirectAbsErr7p95DataSceneIrrMetric.txt";
//const std::string loadFromFileName = "DirectAbsErr4p25DataSceneIrrMetric.txt";
const std::string loadFromFileName = "DirectAbsErr2p1DataSceneIrrMetric.txt";

const char kShaderFile[] = "RenderPasses/AdaptiveSHDemo/AdaptiveGridShader.slang";

#else
//const std::string loadFromFileName = "U32DataScene.txt";
const std::string loadFromFileName = "U32DataScene_4096spp_4spt.txt";
//const std::string loadFromFileName = "U64DataScene.txt";
const char kShaderFile[] = "RenderPasses/AdaptiveSHDemo/UniformGridShader.slang";
#endif

const char kDynamicPassFile[] = "RenderPasses/AdaptiveSHDemo/DynamicPass.slang";
const char kDynamicFilterFile[] = "RenderPasses/AdaptiveSHDemo/DynamicFilter.cs.slang";
const char kCompositePassFile[] = "RenderPasses/AdaptiveSHDemo/CompositeDynamic.slang";
extern "C" FALCOR_API_EXPORT void registerPlugin(Falcor::PluginRegistry& registry)
{
    registry.registerClass<RenderPass, AdaptiveSHDemo>();
}

AdaptiveSHDemo::AdaptiveSHDemo(ref<Device> pDevice, const Properties& props) : RenderPass(pDevice)
{
    mpFbo = Fbo::create(mpDevice);
    Sampler::Desc samplerDesc;
    samplerDesc.setFilterMode(TextureFilteringMode::Linear, TextureFilteringMode::Linear, TextureFilteringMode::Linear);
    mpLinearSampler = mpDevice->createSampler(samplerDesc);
}

Properties AdaptiveSHDemo::getProperties() const
{
    return {};
}

RenderPassReflection AdaptiveSHDemo::reflect(const CompileData& compileData)
{
    // Define the required resources here
    RenderPassReflection reflector;
    const uint2 sz = RenderPassHelpers::calculateIOSize(mOutputSizeSelection, mFixedOutputSize, compileData.defaultTexDims);
    // REMARK MSAA is set via texture sample count. Note that all fbo attachment have to have same sample count.
    reflector.addOutput("output", "Color").texture2D(sz.x, sz.y, 4).format(ResourceFormat::RGBA32Float);
    //reflector.addOutput("output", "Color").texture2D(mLightmapWidth, mLightmapHeight, 4).format(ResourceFormat::RGBA32Float);
    reflector.addOutput("depth", "Depth buffer")
        .format(ResourceFormat::D32Float)
        .bindFlags(ResourceBindFlags::DepthStencil)
        .texture2D(sz.x, sz.y, 4);

    reflector.addInternal("dynamicRaw", "Dynamic raw lighting")
        .format(ResourceFormat::RGBA32Float)
        .bindFlags(ResourceBindFlags::RenderTarget | ResourceBindFlags::ShaderResource)
        .texture2D(sz.x, sz.y);

    reflector.addOutput("dynamicFiltered", "Dynamic filtered lighting")
        .format(ResourceFormat::RGBA32Float)
        .bindFlags(ResourceBindFlags::UnorderedAccess | ResourceBindFlags::ShaderResource)
        .texture2D(sz.x, sz.y);

    reflector.addInternal("dynamicMask", "Dynamic mask")
        .format(ResourceFormat::R8Unorm)
        .bindFlags(ResourceBindFlags::RenderTarget | ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess)
        .texture2D(sz.x, sz.y);

    return reflector;
}

void AdaptiveSHDemo::execute(RenderContext* pRenderContext, const RenderData& renderData)
{
    static uint32_t frameIndex = 0;
    static constexpr uint32_t kWarmupFrames = 100;
    static constexpr uint32_t kMeasureFrames = 1000;

    static auto lastTime = std::chrono::high_resolution_clock::now();
    static std::vector<double> frameTimesMs;

    // STEP 1: Rasterize Scene into UV Space
    auto pTargetFbo = renderData.getTexture("output");
    const float4 clearColor(0, 0, 0, 1);
    mpFbo->attachColorTarget(pTargetFbo, 0);
    // Update frame dimension based on render pass output.
    auto pDepth = renderData.getTexture("depth");
    pRenderContext->clearDsv(pDepth->getDSV().get(), 1.f, 0);
    mpFbo->attachDepthStencilTarget(pDepth);
    pRenderContext->clearFbo(mpFbo.get(), clearColor, 1.0f, 0, FboAttachmentType::Color);
    if (mpScene) {

        // ------------------------------------------------------------------
        // PASS 1: STATIC GEOMETRY (lightmaps)
        // ------------------------------------------------------------------
        auto applyVar = mpStaticVars->getRootVar();
        bindDataSceneData(applyVar);

#if CURRENT_PROBE_MODE == PROBE_MODE_ADAPTIVE
        applyVar["gCornerBuffer"] = mAdaptiveProbeVolume->getCornerBuffer();
        applyVar["gProbeBuffer"] = mAdaptiveProbeVolume->getProbeBuffer();

        applyVar["gSeedProbeIndices"] = mAdaptiveProbeVolume->getSeedProbeIndexBuffer();

        applyVar["SeedGridCB"]["gUseSeedGrid"] = mAdaptiveProbeVolume->getUseSeedGrid() ? 1u : 0u;
        applyVar["SeedGridCB"]["gSeedMinPoint"] = mAdaptiveProbeVolume->getSeedMinPoint();
        applyVar["SeedGridCB"]["gSeedCellSize"] = mAdaptiveProbeVolume->getSeedCellSize();
        applyVar["SeedGridCB"]["gSeedResolution"] = mAdaptiveProbeVolume->getSeedResolution();

#else
        mUniformProbeVolume->bindShaderData(applyVar);
#endif

        mpScene->rasterize(pRenderContext, mpGraphicsState.get(), mpStaticVars.get(), mpRasterState, mpRasterState);

        auto now = std::chrono::high_resolution_clock::now();
        double frameMs = std::chrono::duration<double, std::milli>(now - lastTime).count();
        lastTime = now;

        if (frameIndex >= kWarmupFrames && frameIndex < kWarmupFrames + kMeasureFrames)
        {
            frameTimesMs.push_back(frameMs);
        }

        if (frameIndex == kWarmupFrames + kMeasureFrames)
        {
            double sum = 0.0;
            for (double v : frameTimesMs) sum += v;

            double meanMs = sum / double(frameTimesMs.size());
            double meanFps = 1000.0 / meanMs;

            double variance = 0.0;
            for (double v : frameTimesMs)
            {
                double d = v - meanMs;
                variance += d * d;
            }
            variance /= double(frameTimesMs.size());
            double stdMs = std::sqrt(variance);
            uint64_t probeCount = 0;
#if CURRENT_PROBE_MODE == PROBE_MODE_ADAPTIVE
            std::string modeName = "Adaptive";
            probeCount = mAdaptiveProbeVolume->getProbes().size();
#else
            std::string modeName = "Uniform";
            const uint3 res = mUniformProbeVolume->getCellResolution();
            probeCount = uint64_t(res.x) * uint64_t(res.y) * uint64_t(res.z);
#endif

            //std::ofstream out("AdaptiveSHDemo_RuntimeFPS.csv", std::ios::app);
            //out << modeName << ","
            //    << loadFromFileName << ","
            //    << probeCount << ","
            //    << kWarmupFrames << ","
            //    << kMeasureFrames << ","
            //    << std::fixed << std::setprecision(4)
            //    << meanMs << ","
            //    << stdMs << ","
            //    << meanFps << "\n";

            //out.close();

            //logInfo("Runtime measurement finished: file={}, mean {:.4f} ms, std {:.4f} ms, FPS {:.2f}",
            //    loadFromFileName, meanMs, stdMs, meanFps);
        }

        frameIndex++;

        if (mbShowAdaptiveGrid)
        {
            mpProbeVisualizePass->setCameraData(
                mpScene->getCamera()->getViewProjMatrix()
            );

            mpProbeVisualizePass->execute(pRenderContext, mpFbo);
        }
    }
}

void AdaptiveSHDemo::renderUI(Gui::Widgets& widget) {
    widget.text("Loaded probe file: " + loadFromFileName);
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

void AdaptiveSHDemo::loadLightmaps()
{
    // Load as a 2D texture. Falcor handles EXR (HDR) automatically.
    // We set loadAsSrgb to false because lightmaps contain linear radiance data.
    //auto loadOne = [&](size_t idx, ref<Texture>& dst, const std::string& debugName)
    //    {
    //        if (idx >= mBakeTargets.size())
    //        {
    //            logWarning("Bake target index {} is out of range.", idx);
    //            return;
    //        }

    //        const auto& target = mBakeTargets[idx];

    //        dst = Texture::createFromFile(
    //            mpDevice,
    //            target.outputPath,
    //            true,
    //            false,
    //            ResourceBindFlags::ShaderResource
    //        );

    //        if (dst)
    //        {
    //            dst->setName(debugName);
    //            logInfo("Successfully loaded lightmap '{}' from {}", target.name, target.outputPath);
    //        }
    //        else
    //        {
    //            logWarning("Failed to load lightmap '{}' from {}", target.name, target.outputPath);
    //        }
    //    };

    //loadOne(0, mpFloorLightmap, "FloorLightmap");
    //loadOne(1, mpLeftWallLightmap, "LeftWallLightmap");
    //loadOne(2, mpRightWallLightmap, "RightWallLightmap");
    //loadOne(3, mpRoofLeftLightmap, "RoofLeftLightmap");
    //loadOne(4, mpRoofRightLightmap, "RoofRightLightmap");

    //loadOne(5, mpPillar0Lightmap, "Pillar0Lightmap");
    //loadOne(6, mpPillar1Lightmap, "Pillar1Lightmap");
    //loadOne(7, mpPillar2Lightmap, "Pillar2Lightmap");
    //loadOne(8, mpPillar3Lightmap, "Pillar3Lightmap");
    //loadOne(9, mpPillar4Lightmap, "Pillar4Lightmap");
    //loadOne(10, mpPillar5Lightmap, "Pillar5Lightmap");
    //loadOne(11, mpPillar6Lightmap, "Pillar6Lightmap");
    //loadOne(12, mpPillar7Lightmap, "Pillar7Lightmap");
}

void AdaptiveSHDemo::setScene(RenderContext* pRenderContext, const ref<Scene>& pScene)
{
    mpScene = pScene;
    if (mpScene)
    {
        std::ifstream check("AdaptiveSHDemo_RuntimeFPS.csv");
        if (!check.good())
        {
            std::ofstream out("AdaptiveSHDemo_RuntimeFPS.csv");
            out << "Mode,ProbeFile,ProbeCount,WarmupFrames,MeasuredFrames,MeanFrameMs,StdFrameMs,MeanFPS\n";
        }

        setupDataSceneBakeTargets();

        //init probe visual pass
        mpProbeVisualizePass = ProbeVisualizePass::create(mpDevice, mpScene->getSceneDefines());
#if CURRENT_PROBE_MODE == PROBE_MODE_ADAPTIVE
        mAdaptiveProbeVolume = AdaptiveProbeVolume::create(mpDevice);

            mAdaptiveProbeVolume->loadFromFile(loadFromFileName);
            mAdaptiveProbeVolume->uploadToGPU();
            mpProbeVisualizePass->setVolumeData(mAdaptiveProbeVolume->getProbes());
#else
        mUniformProbeVolume = UniformProbeVolume::create(mpDevice);
            mUniformProbeVolume->loadFromFile(loadFromFileName);
            mUniformProbeVolume->uploadToGPU();
            mpProbeVisualizePass->setUniformVolumeData(
                mUniformProbeVolume->getMinPoint(),
                mUniformProbeVolume->getCellSize(),
                mUniformProbeVolume->getCellResolution()
            );
#endif
      

        ProgramDesc previewDesc;
        previewDesc.addShaderModules(mpScene->getShaderModules());
        previewDesc.addShaderLibrary(kShaderFile)
            .vsEntry("vsMain")
            .psEntry("psMain");

        mpStaticProgram = Program::create(mpDevice, previewDesc, mpScene->getSceneDefines());
        mpStaticVars = ProgramVars::create(mpDevice, mpStaticProgram->getReflector());

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
        mpGraphicsState->setProgram(mpStaticProgram);
        mpGraphicsState->setRasterizerState(mpRasterState);
        mpGraphicsState->setFbo(mpFbo);
        mpGraphicsState->setDepthStencilState(pDsState);

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

        ProgramDesc dynamicDesc;
        dynamicDesc.addShaderModules(mpScene->getShaderModules());
        dynamicDesc.addShaderLibrary(kDynamicPassFile)
            .vsEntry("vsMain")
            .psEntry("psMain");
        dynamicDesc.addTypeConformances(mpScene->getTypeConformances());

        DefineList dynamicDefines = mpScene->getSceneDefines();
        DefineList lightRelatedDefines;
        lightRelatedDefines.add("USE_ANALYTIC_LIGHTS", mpScene->useAnalyticLights() ? "1" : "0");
        lightRelatedDefines.add("USE_EMISSIVE_LIGHTS", mpScene->useEmissiveLights() ? "1" : "0");
        dynamicDefines.add(lightRelatedDefines);

        if (mpEmissiveSampler)
        {
            dynamicDefines.add(mpEmissiveSampler->getDefines());
        }

        mpDynamicProgram = Program::create(mpDevice, dynamicDesc, dynamicDefines);
        mpDynamicVars = ProgramVars::create(mpDevice, mpDynamicProgram->getReflector());

        mpFilterPass = ComputePass::create(mpDevice, kDynamicFilterFile, "main");

        ProgramDesc compositeDesc;
        compositeDesc.addShaderLibrary(kCompositePassFile)
            .vsEntry("vsMain")
            .psEntry("psMain");

        mpCompositeProgram = Program::create(mpDevice, compositeDesc);
        mpCompositeVars = ProgramVars::create(mpDevice, mpCompositeProgram->getReflector());

        BlendState::Desc blendDesc;
        blendDesc.setRtBlend(0, true)
            .setRtParams(0,
                BlendState::BlendOp::Add,
                BlendState::BlendOp::Add,
                BlendState::BlendFunc::SrcAlpha,
                BlendState::BlendFunc::OneMinusSrcAlpha,
                BlendState::BlendFunc::One,
                BlendState::BlendFunc::OneMinusSrcAlpha);

        DepthStencilState::Desc dsComposite;
        dsComposite.setDepthEnabled(false);

        mpCompositeState = GraphicsState::create(mpDevice);
        mpCompositeState->setProgram(mpCompositeProgram);
        mpCompositeState->setBlendState(BlendState::create(blendDesc));
        mpCompositeState->setDepthStencilState(DepthStencilState::create(dsComposite));
        mpCompositeState->setRasterizerState(mpRasterState);

        mpEmptyVao = Vao::create(Vao::Topology::TriangleStrip);
        mpCompositeState->setVao(mpEmptyVao);

        //loadLightmaps();
        loadDataSceneLightmaps();
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

void AdaptiveSHDemo::loadDataSceneLightmaps()
{
    auto loadOne = [&](size_t idx, ref<Texture>& dst, const std::string& debugName)
        {
            if (idx >= mBakeTargets.size())
            {
                logWarning("Bake target index {} is out of range.", idx);
                return;
            }

            const auto& target = mBakeTargets[idx];

            dst = Texture::createFromFile(
                mpDevice,
                target.outputPath,
                true,
                false,
                ResourceBindFlags::ShaderResource
            );

            if (dst)
            {
                dst->setName(debugName);
                logInfo("Successfully loaded data-scene lightmap '{}' from {}", target.name, target.outputPath);
            }
            else
            {
                logWarning("Failed to load data-scene lightmap '{}' from {}", target.name, target.outputPath);
            }
        };

    loadOne(0, mpDataFloorLightmap, "DataFloorLightmap");
    loadOne(1, mpDataCeilingLightmap, "DataCeilingLightmap");
    loadOne(2, mpDataLeftWallLightmap, "DataLeftWallLightmap");
    loadOne(3, mpDataRightWallLightmap, "DataRightWallLightmap");
    loadOne(4, mpDataBackWallLightmap, "DataBackWallLightmap");
    loadOne(5, mpDataFrontWallLightmap, "DataFrontWallLightmap");

    loadOne(6, mpDataTallBoxALightmap, "DataTallBoxALightmap");
    loadOne(7, mpDataWideBoxBLightmap, "DataWideBoxBLightmap");
    loadOne(8, mpDataBlockCLightmap, "DataBlockCLightmap");
    loadOne(9, mpDataLowBoxDLightmap, "DataLowBoxDLightmap");
    loadOne(10, mpDataThinSlabELightmap, "DataThinSlabELightmap");
    loadOne(11, mpDataThinSlabFLightmap, "DataThinSlabFLightmap");
    loadOne(12, mpDataTallBoxJLightmap, "DataTallBoxJLightmap");
    loadOne(13, mpDataShortBoxLLightmap, "DataShortBoxLLightmap");
}

void AdaptiveSHDemo::bindDataSceneData(ShaderVar applyVar)
{
    applyVar["gLinearSampler"] = mpLinearSampler;

    applyVar["gFloorLightmap"] = mpDataFloorLightmap;
    applyVar["gCeilingLightmap"] = mpDataCeilingLightmap;
    applyVar["gLeftWallLightmap"] = mpDataLeftWallLightmap;
    applyVar["gRightWallLightmap"] = mpDataRightWallLightmap;
    applyVar["gBackWallLightmap"] = mpDataBackWallLightmap;
    applyVar["gFrontWallLightmap"] = mpDataFrontWallLightmap;

    applyVar["gTallBoxALightmap"] = mpDataTallBoxALightmap;
    applyVar["gWideBoxBLightmap"] = mpDataWideBoxBLightmap;
    applyVar["gBlockCLightmap"] = mpDataBlockCLightmap;
    applyVar["gLowBoxDLightmap"] = mpDataLowBoxDLightmap;
    applyVar["gThinSlabELightmap"] = mpDataThinSlabELightmap;
    applyVar["gThinSlabFLightmap"] = mpDataThinSlabFLightmap;
    applyVar["gTallBoxJLightmap"] = mpDataTallBoxJLightmap;
    applyVar["gShortBoxLLightmap"] = mpDataShortBoxLLightmap;

    applyVar["PerFrameCB"]["gFloorInstanceID"] = 0;
    applyVar["PerFrameCB"]["gCeilingInstanceID"] = 1;
    applyVar["PerFrameCB"]["gLeftWallInstanceID"] = 2;
    applyVar["PerFrameCB"]["gRightWallInstanceID"] = 3;
    applyVar["PerFrameCB"]["gBackWallInstanceID"] = 4;
    applyVar["PerFrameCB"]["gFrontWallInstanceID"] = 5;

    applyVar["PerFrameCB"]["gTallBoxAInstanceID"] = 6;
    applyVar["PerFrameCB"]["gWideBoxBInstanceID"] = 7;
    applyVar["PerFrameCB"]["gBlockCInstanceID"] = 8;
    applyVar["PerFrameCB"]["gLowBoxDInstanceID"] = 9;
    applyVar["PerFrameCB"]["gThinSlabEInstanceID"] = 10;
    applyVar["PerFrameCB"]["gThinSlabFInstanceID"] = 11;
    applyVar["PerFrameCB"]["gTallBoxJInstanceID"] = 12;
    applyVar["PerFrameCB"]["gShortBoxLInstanceID"] = 13;

    applyVar["PerFrameCB"]["gTallBoxACenterW"] = float3(-3.8f, 1.3f, -2.8f);
    applyVar["PerFrameCB"]["gWideBoxBCenterW"] = float3(-1.9f, 0.5f, 2.2f);
    applyVar["PerFrameCB"]["gBlockCCenterW"] = float3(2.8f, 0.8f, -2.6f);
    applyVar["PerFrameCB"]["gLowBoxDCenterW"] = float3(1.4f, 0.35f, 2.8f);
    applyVar["PerFrameCB"]["gThinSlabECenterW"] = float3(-0.4f, 1.5f, -0.8f);
    applyVar["PerFrameCB"]["gThinSlabFCenterW"] = float3(3.4f, 1.2f, 0.8f);
    applyVar["PerFrameCB"]["gTallBoxJCenterW"] = float3(-4.2f, 1.5f, 2.8f);
    applyVar["PerFrameCB"]["gShortBoxLCenterW"] = float3(-0.8f, 0.45f, -3.3f);

    applyVar["PerFrameCB"]["gTallBoxAHalfExtentW"] = float3(0.45f, 1.30f, 0.45f);
    applyVar["PerFrameCB"]["gWideBoxBHalfExtentW"] = float3(0.90f, 0.50f, 0.60f);
    applyVar["PerFrameCB"]["gBlockCHalfExtentW"] = float3(0.72f, 0.80f, 0.72f);
    applyVar["PerFrameCB"]["gLowBoxDHalfExtentW"] = float3(1.05f, 0.35f, 1.00f);
    applyVar["PerFrameCB"]["gThinSlabEHalfExtentW"] = float3(0.175f, 1.50f, 0.70f);
    applyVar["PerFrameCB"]["gThinSlabFHalfExtentW"] = float3(0.36f, 1.20f, 0.85f);
    applyVar["PerFrameCB"]["gTallBoxJHalfExtentW"] = float3(0.75f, 1.25f, 0.75f);
    applyVar["PerFrameCB"]["gShortBoxLHalfExtentW"] = float3(0.75f, 0.45f, 0.75f);
}

void AdaptiveSHDemo::setupDataSceneBakeTargets()
{
    mBakeTargets =
    {
        // Room shell.
        //{ "Floor",     0, 1024, 1024, "BakedData_Floor.exr"     },
        { "Floor",     0, 2048, 2048, "BakedData_Floor.exr"     },
        { "Ceiling",   1, 1024, 1024, "BakedData_Ceiling.exr"   },
        { "LeftWall",  2, 1024,  512, "BakedData_LeftWall.exr"  },
        { "RightWall", 3, 1024,  512, "BakedData_RightWall.exr" },
        { "BackWall",  4, 1024,  512, "BakedData_BackWall.exr"  },
        { "FrontWall", 5, 1024,  512, "BakedData_FrontWall.exr" },

        //Box targets.
        //halfExtent = scaling * 0.5
        //rotationEulerDeg = same as SceneForData.pyscene

       { "TallBoxA", 6, 512, 512, "BakedData_TallBoxA.exr",
           BakeTargetType::Pillar,
           float3(-3.8f, 1.3f, -2.8f),
           float3(0.45f, 1.30f, 0.45f),
           float3(0.0f, 0.0f, 0.0f) },

       { "WideBoxB", 7, 1024, 1024, "BakedData_WideBoxB.exr",
           BakeTargetType::Pillar,
           float3(-1.9f, 0.5f, 2.2f),
           float3(0.90f, 0.50f, 0.60f),
           float3(0.0f, 0.0f, 0.0f) },

       { "BlockC", 8, 512, 512, "BakedData_BlockC.exr",
           BakeTargetType::Pillar,
           float3(2.8f, 0.8f, -2.6f),
           float3(0.55f, 0.80f, 0.55f),
           float3(0.0f, 18.0f, 0.0f) },

       { "LowBoxD", 9, 512, 512, "BakedData_LowBoxD.exr",
           BakeTargetType::Pillar,
           float3(1.4f, 0.35f, 2.8f),
           float3(0.80f, 0.35f, 0.70f),
           float3(0.0f, -22.0f, 0.0f) },

       { "ThinSlabE", 10, 512, 512, "BakedData_ThinSlabE.exr",
           BakeTargetType::Pillar,
           float3(-0.4f, 1.5f, -0.8f),
           float3(0.175f, 1.50f, 0.70f),
           float3(0.0f, 0.0f, 0.0f) },

       { "ThinSlabF", 11, 512, 512, "BakedData_ThinSlabF.exr",
           BakeTargetType::Pillar,
           float3(3.4f, 1.2f, 0.8f),
           float3(0.175f, 1.20f, 0.80f),
           float3(0.0f, 12.0f, 0.0f) },

       { "TallBoxJ", 12, 512, 512, "BakedData_TallBoxJ.exr",
           BakeTargetType::Pillar,
           float3(-4.2f, 1.5f, 2.8f),
           float3(0.35f, 1.10f, 0.35f),
           float3(0.0f, -18.0f, 8.0f) },

       { "ShortBoxL", 13, 512, 512, "BakedData_ShortBoxL.exr",
           BakeTargetType::Pillar,
           float3(-0.8f, 0.45f, -3.3f),
           float3(0.50f, 0.45f, 0.50f),
           float3(0.0f, 40.0f, 0.0f) },
    };
}
