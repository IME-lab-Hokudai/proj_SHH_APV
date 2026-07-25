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
#include "BakeLightMapSingle.h"
#include "Rendering/Lights/EmissivePowerSampler.h"
#include "Rendering/Lights/EmissiveUniformSampler.h"
#include "BakeDataStructures.slang"
#include "Utils/Logger.h"
const char kBakingFile[] = "RenderPasses/BakeLightMapSingle/LightmapBakingSingle.rt.slang";
const char kUVrasterFile[] = "RenderPasses/BakeLightMapSingle/UVPassSingle.slang";
const char kUVrasterPillarFile[] = "RenderPasses/BakeLightMapSingle/UVPassPillar.slang";
const char kExtractFile[] = "RenderPasses/BakeLightMapSingle/ExtractTexelsSingle.cs.slang";
const char kNormalizeFile[] = "RenderPasses/BakeLightMapSingle/NormalizeLightmapSingle.cs.slang";
const char kShaderFile[] = "RenderPasses/BakeLightMapSingle/ApplyLightmapSingle.slang";
const char kDataShaderFile[] = "RenderPasses/BakeLightMapSingle/ApplyLightmapDataScene.slang";

extern "C" FALCOR_API_EXPORT void registerPlugin(Falcor::PluginRegistry& registry)
{
    registry.registerClass<RenderPass, BakeLightMapSingle>();
}

BakeLightMapSingle::BakeLightMapSingle(ref<Device> pDevice, const Properties& props) : RenderPass(pDevice)
{
    mpFbo = Fbo::create(mpDevice);
    Sampler::Desc samplerDesc;
    samplerDesc.setFilterMode(TextureFilteringMode::Linear, TextureFilteringMode::Linear, TextureFilteringMode::Linear);
    mpLinearSampler = mpDevice->createSampler(samplerDesc);
}

Properties BakeLightMapSingle::getProperties() const
{
    return {};
}

RenderPassReflection BakeLightMapSingle::reflect(const CompileData& compileData)
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
    return reflector;
}

void BakeLightMapSingle::execute(RenderContext* pRenderContext, const RenderData& renderData)
{
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
        if (mbloadLightMap) {
            auto applyVar = mpVars->getRootVar();
            //bindDataSceneData(applyVar);
            bindCornellData(applyVar);
            mpScene->rasterize(pRenderContext, mpGraphicsState.get(), mpVars.get(), mpRasterState, mpRasterState);
        }
        else {
            if (mCurrentTargetIndex < mBakeTargets.size())
            {
                const auto& target = mBakeTargets[mCurrentTargetIndex];

                if (target.type == BakeTargetType::Pillar)
                {
                    mpGraphicsState->setProgram(mpProgramPillar);
                    mpVars = ProgramVars::create(mpDevice, mpProgramPillar->getReflector());
                }
                else
                {
                    mpGraphicsState->setProgram(mpProgram);
                    mpVars = ProgramVars::create(mpDevice, mpProgram->getReflector());
                }
                // -------------------------------------------------------------------------
                // PREPARE CURRENT TARGET (only once per target)
                // -------------------------------------------------------------------------
                if (mNeedsPreparation)
                {
                    mLightmapWidth = target.width;
                    mLightmapHeight = target.height;
                    GraphicsState::Viewport vp(0.0f, 0.0f, (float)mLightmapWidth, (float)mLightmapHeight, 0.0f, 1.0f);
                    mpGraphicsState->setViewport(0, vp);
                    ref<Texture> pPosTex = mpDevice->createTexture2D(
                        mLightmapWidth, mLightmapHeight, ResourceFormat::RGBA32Float,
                        1, 1, nullptr,
                        ResourceBindFlags::RenderTarget | ResourceBindFlags::ShaderResource
                    );
                    pPosTex->setName("PosTex");

                    ref<Texture> pNormTex = mpDevice->createTexture2D(
                        mLightmapWidth, mLightmapHeight, ResourceFormat::RGBA32Float,
                        1, 1, nullptr,
                        ResourceBindFlags::RenderTarget | ResourceBindFlags::ShaderResource
                    );
                    pNormTex->setName("NormTex");

                    mpUVFbo = Fbo::create(mpDevice);
                    mpUVFbo->attachColorTarget(pPosTex, 0);
                    mpUVFbo->attachColorTarget(pNormTex, 1);

                    uint32_t totalTexels = mLightmapWidth * mLightmapHeight;

                    mpTexelBuffer = mpDevice->createStructuredBuffer(sizeof(TexelSample), totalTexels);
                    mpTexelBuffer->setName("TexelBuffer");

                    mpCounterBuffer = mpDevice->createBuffer(
                        sizeof(uint32_t),
                        ResourceBindFlags::UnorderedAccess,
                        MemoryType::DeviceLocal
                    );
                    mpCounterBuffer->setName("CounterBuffer");

                    mpAccumBuffer = mpDevice->createStructuredBuffer(sizeof(float4), totalTexels);
                    mpAccumBuffer->setName("AccumBuffer");

                    mpResultTex = mpDevice->createTexture2D(
                        mLightmapWidth, mLightmapHeight, ResourceFormat::RGBA16Float,
                        1, 1, nullptr,
                        ResourceBindFlags::ShaderResource | ResourceBindFlags::RenderTarget | ResourceBindFlags::UnorderedAccess
                    );
                    mpResultTex->setName("BakeResultTex");

                    mpGraphicsState->setFbo(mpUVFbo);
                    pRenderContext->clearFbo(mpUVFbo.get(), float4(0), 1.0f, 0, FboAttachmentType::Color);

                    auto var = mpVars->getRootVar();
                    var["BakeCB"]["gReceiverInstanceID"] = target.instanceID;
                    var["BakeCB"]["gPillarCenterW"] = target.pillarCenterW;
                    var["BakeCB"]["gPillarHalfExtentW"] = target.pillarHalfExtentW;
                    var["BakeCB"]["gPillarRotationEulerDeg"] = target.pillarRotationEulerDeg;
                    mpScene->rasterize(pRenderContext, mpGraphicsState.get(), mpVars.get(), mpRasterState, mpRasterState);

                    pRenderContext->clearUAV(mpCounterBuffer->getUAV().get(), uint4(0));

                    auto extractVar = mpExtractPass->getRootVar();
                    extractVar["gPosW"] = mpUVFbo->getColorTexture(0);
                    extractVar["gNormW"] = mpUVFbo->getColorTexture(1);
                    extractVar["gTexelSamples"] = mpTexelBuffer;
                    extractVar["gCounter"] = mpCounterBuffer;
                    extractVar["CB"]["gResolution"] = uint2(mLightmapWidth, mLightmapHeight);

                    mpExtractPass->execute(pRenderContext, mLightmapWidth, mLightmapHeight);

                    uint32_t extractedCount = 0;
                    mpCounterBuffer->getBlob(&extractedCount, 0, sizeof(uint32_t));
                    mNumExtractedTexels = extractedCount;

                    logInfo("Prepared target '{}' (instanceID={}), extracted {} texels.",
                        target.name, target.instanceID, mNumExtractedTexels);

                    mCurrentSample = 0;
                    mNeedsPreparation = false;

                    if (mNumExtractedTexels == 0)
                    {
                        logWarning("Target '{}' has no extracted texels. Skipping.", target.name);
                        mCurrentTargetIndex++;
                        //mCurrentTargetIndex = 0;
                        mNeedsPreparation = true;
                        return;
                    }
                }

                // -------------------------------------------------------------------------
                // ITERATIVE BAKING (Ray Tracing Accumulation)
                // -------------------------------------------------------------------------
                if (mCurrentSample < mTotalSamples)
                {
                    if (mCurrentSample == 0)
                    {
                        pRenderContext->clearUAV(mpAccumBuffer->getUAV().get(), float4(0));
                    }

                    auto rtVar = mpRtVars->getRootVar();
                    rtVar["gTexelSamples"] = mpTexelBuffer;
                    rtVar["gIrradianceAccum"] = mpAccumBuffer;
                    rtVar["PerFrameCB"]["sampleIndex"] = mCurrentSample;
                    rtVar["PerFrameCB"]["numTexels"] = mNumExtractedTexels;
                    rtVar["PerFrameCB"]["bias"] = 0.01f;

                    if (mpEmissiveSampler)
                    {
                        mpEmissiveSampler->bindShaderData(rtVar["PerFrameCB"]["emissiveSampler"]);
                    }

                    uint32_t threadCount = mNumExtractedTexels;
                    mpScene->raytrace(pRenderContext, mpRtProgram.get(), mpRtVars, uint3(threadCount, 1, 1));

                    mCurrentSample++;

                    // ---------------------------------------------------------------------
                    // FINALIZATION
                    // ---------------------------------------------------------------------
                    if (mCurrentSample == mTotalSamples)
                    {
                        auto normalizeVar = mpNormalizePass->getRootVar();
                        normalizeVar["gAccumBuffer"] = mpAccumBuffer;
                        normalizeVar["gOutput"] = mpResultTex;
                        normalizeVar["CB"]["gTotalSamples"] = mTotalSamples;
                        normalizeVar["CB"]["gWidth"] = mLightmapWidth;

                        mpNormalizePass->execute(pRenderContext, mLightmapWidth, mLightmapHeight);

                        if (mpResultTex)
                        {
                            mpResultTex->captureToFile(
                                0,
                                0,
                                target.outputPath,
                                Bitmap::FileFormat::ExrFile,
                                Bitmap::ExportFlags::Uncompressed
                            );
                            logInfo("Lightmap saved to: {}", target.outputPath);
                        }

                        logInfo("Finished baking '{}' ({} samples).", target.name, mTotalSamples);

                        mCurrentTargetIndex++;
                        mCurrentSample = 0;
                        mNumExtractedTexels = 0;
                        mNeedsPreparation = true;
                    }
                }
            }
        }
    }
}

void BakeLightMapSingle::renderUI(Gui::Widgets& widget) {}

void BakeLightMapSingle::loadLightmaps()
{
    // Load as a 2D texture. Falcor handles EXR (HDR) automatically.
    // We set loadAsSrgb to false because lightmaps contain linear radiance data.
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
                logInfo("Successfully loaded lightmap '{}' from {}", target.name, target.outputPath);
            }
            else
            {
                logWarning("Failed to load lightmap '{}' from {}", target.name, target.outputPath);
            }
        };

    loadOne(0, mpFloorLightmap, "FloorLightmap");
    loadOne(1, mpLeftWallLightmap, "LeftWallLightmap");
    loadOne(2, mpRightWallLightmap, "RightWallLightmap");
    loadOne(3, mpRoofLeftLightmap, "RoofLeftLightmap");
    loadOne(4, mpRoofRightLightmap, "RoofRightLightmap");

    loadOne(5, mpPillar0Lightmap, "Pillar0Lightmap");
    loadOne(6, mpPillar1Lightmap, "Pillar1Lightmap");
    loadOne(7, mpPillar2Lightmap, "Pillar2Lightmap");
    loadOne(8, mpPillar3Lightmap, "Pillar3Lightmap");
    loadOne(9, mpPillar4Lightmap, "Pillar4Lightmap");
    loadOne(10, mpPillar5Lightmap, "Pillar5Lightmap");
    loadOne(11, mpPillar6Lightmap, "Pillar6Lightmap");
    loadOne(12, mpPillar7Lightmap, "Pillar7Lightmap");
}

void BakeLightMapSingle::setScene(RenderContext* pRenderContext, const ref<Scene>& pScene)
{
    mpScene = pScene;
    if (mpScene)
    {
        //setupDataSceneBakeTargets();
        setupCornellBakeTargets();
        if (mbloadLightMap) {
            ProgramDesc previewDesc;
            previewDesc.addShaderModules(mpScene->getShaderModules());
            //previewDesc.addShaderLibrary(kShaderFile)
            //    .vsEntry("vsMain")
            //    .psEntry("psMain");
            previewDesc.addShaderLibrary(kDataShaderFile)
                .vsEntry("vsMain")
                .psEntry("psMain");

            mpProgram = Program::create(mpDevice, previewDesc, mpScene->getSceneDefines());
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

            //loadLightmaps();
            //loadDataSceneLightmaps();
            loadCornellLightmaps();
        }
        else {
            //logInfo("Geometry instance count = {}", mpScene->getGeometryInstanceCount());

            //for (uint32_t instanceID = 0; instanceID < mpScene->getGeometryInstanceCount(); ++instanceID)
            //{
            //    const auto& inst = mpScene->getGeometryInstance(instanceID);
            //    logInfo("instanceID={} meshID={} materialID={} globalMatrixID={}",
            //        (uint32_t)instanceID,
            //        (uint32_t)inst.geometryID,
            //        (uint32_t)inst.materialID,
            //        (uint32_t)inst.globalMatrixID);
            //}

            ProgramDesc desc;
            desc.addShaderModules(mpScene->getShaderModules());
            desc.addShaderLibrary(kUVrasterFile)
                .vsEntry("vsMain")  // Vertex shader entry point
                .psEntry("psMain"); // Pixel shader entry point;
            mpProgram = Program::create(mpDevice, desc, mpScene->getSceneDefines());

            ProgramDesc descPillar;
            descPillar.addShaderModules(mpScene->getShaderModules());
            descPillar.addShaderLibrary(kUVrasterPillarFile)
                .vsEntry("vsMain")  // Vertex shader entry point
                .psEntry("psMain"); // Pixel shader entry point;
            mpProgramPillar = Program::create(mpDevice, descPillar, mpScene->getSceneDefines());

            // rasterizer state
            RasterizerState::Desc rasterDesc;
            rasterDesc.setFillMode(RasterizerState::FillMode::Solid);
            rasterDesc.setCullMode(RasterizerState::CullMode::None);

            mpRasterState = RasterizerState::create(rasterDesc);

            // default depth stencil state
            DepthStencilState::Desc dsDesc;
            dsDesc.setDepthEnabled(false); // Add this line

            ref<DepthStencilState> pDsState = DepthStencilState::create(dsDesc);

            mpGraphicsState = GraphicsState::create(mpDevice);
            mpGraphicsState->setRasterizerState(mpRasterState);
            mpGraphicsState->setDepthStencilState(pDsState);

            mpExtractPass = ComputePass::create(mpDevice, kExtractFile, "main", mpScene->getSceneDefines());
            mpNormalizePass = ComputePass::create(mpDevice, kNormalizeFile, "main", mpScene->getSceneDefines());
            ProgramDesc rtProgDesc;
            rtProgDesc.addShaderModules(mpScene->getShaderModules());
            rtProgDesc.addShaderLibrary(kBakingFile);
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
}

void BakeLightMapSingle::setupDataSceneBakeTargets()
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

void BakeLightMapSingle::bindDataCorridor(ShaderVar applyVar)
{
    applyVar["gLinearSampler"] = mpLinearSampler; // Standard linear sampler
    applyVar["gFloorLightmap"] = mpFloorLightmap;
    applyVar["gLeftWallLightmap"] = mpLeftWallLightmap;
    applyVar["gRightWallLightmap"] = mpRightWallLightmap;
    applyVar["gRoofLeftLightmap"] = mpRoofLeftLightmap;
    applyVar["gRoofRightLightmap"] = mpRoofRightLightmap;
    applyVar["gPillar0Lightmap"] = mpPillar0Lightmap;
    applyVar["gPillar1Lightmap"] = mpPillar1Lightmap;
    applyVar["gPillar2Lightmap"] = mpPillar2Lightmap;
    applyVar["gPillar3Lightmap"] = mpPillar3Lightmap;
    applyVar["gPillar4Lightmap"] = mpPillar4Lightmap;
    applyVar["gPillar5Lightmap"] = mpPillar5Lightmap;
    applyVar["gPillar6Lightmap"] = mpPillar6Lightmap;
    applyVar["gPillar7Lightmap"] = mpPillar7Lightmap;

    applyVar["PerFrameCB"]["gFloorInstanceID"] = mFloorInstanceID;
    applyVar["PerFrameCB"]["gLeftWallInstanceID"] = mLeftWallInstanceID;
    applyVar["PerFrameCB"]["gRightWallInstanceID"] = mRightWallInstanceID;
    applyVar["PerFrameCB"]["gRoofLeftInstanceID"] = mRoofLeftInstanceID;
    applyVar["PerFrameCB"]["gRoofRightInstanceID"] = mRoofRightInstanceID;
    applyVar["PerFrameCB"]["gPillar0InstanceID"] = mPillar0InstanceID;
    applyVar["PerFrameCB"]["gPillar1InstanceID"] = mPillar1InstanceID;
    applyVar["PerFrameCB"]["gPillar2InstanceID"] = mPillar2InstanceID;
    applyVar["PerFrameCB"]["gPillar3InstanceID"] = mPillar3InstanceID;
    applyVar["PerFrameCB"]["gPillar4InstanceID"] = mPillar4InstanceID;
    applyVar["PerFrameCB"]["gPillar5InstanceID"] = mPillar5InstanceID;
    applyVar["PerFrameCB"]["gPillar6InstanceID"] = mPillar6InstanceID;
    applyVar["PerFrameCB"]["gPillar7InstanceID"] = mPillar7InstanceID;

    applyVar["PerFrameCB"]["gPillar0CenterW"] = float3(7.5f, 5.0f, 5.0f);
    applyVar["PerFrameCB"]["gPillar1CenterW"] = float3(-7.5f, 5.0f, 5.0f);
    applyVar["PerFrameCB"]["gPillar2CenterW"] = float3(7.5f, 5.0f, 15.0f);
    applyVar["PerFrameCB"]["gPillar3CenterW"] = float3(-7.5f, 5.0f, 15.0f);
    applyVar["PerFrameCB"]["gPillar4CenterW"] = float3(7.5f, 5.0f, 25.0f);
    applyVar["PerFrameCB"]["gPillar5CenterW"] = float3(-7.5f, 5.0f, 25.0f);
    applyVar["PerFrameCB"]["gPillar6CenterW"] = float3(7.5f, 5.0f, 35.0f);
    applyVar["PerFrameCB"]["gPillar7CenterW"] = float3(-7.5f, 5.0f, 35.0f);

    applyVar["PerFrameCB"]["gPillar0HalfExtentW"] = float3(0.75f, 5.0f, 0.75f);
    applyVar["PerFrameCB"]["gPillar1HalfExtentW"] = float3(0.75f, 5.0f, 0.75f);
    applyVar["PerFrameCB"]["gPillar2HalfExtentW"] = float3(0.75f, 5.0f, 0.75f);
    applyVar["PerFrameCB"]["gPillar3HalfExtentW"] = float3(0.75f, 5.0f, 0.75f);
    applyVar["PerFrameCB"]["gPillar4HalfExtentW"] = float3(0.75f, 5.0f, 0.75f);
    applyVar["PerFrameCB"]["gPillar5HalfExtentW"] = float3(0.75f, 5.0f, 0.75f);
    applyVar["PerFrameCB"]["gPillar6HalfExtentW"] = float3(0.75f, 5.0f, 0.75f);
    applyVar["PerFrameCB"]["gPillar7HalfExtentW"] = float3(0.75f, 5.0f, 0.75f);
}

void BakeLightMapSingle::bindDataSceneData(ShaderVar applyVar)
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

void BakeLightMapSingle::setupCorridorBakeTargets()
{
    mBakeTargets =
    {
        { "Floor",     0, 1024, 1024, "BakedFloor.exr"     },
        { "LeftWall",  1, 1024,  512, "BakedLeftWall.exr"  },
        { "RightWall", 2, 1024,  512, "BakedRightWall.exr" },
        { "RoofLeft", 11, 1024,  512, "BakedRoofLeft.exr"  },
        { "RoofRight",12, 1024,  512, "BakedRoofRight.exr" },

        { "Pillar0",  3, 512, 512, "BakedPillar0.exr",
            BakeTargetType::Pillar, float3(7.5f, 5.0f,  5.0f), float3(0.75f, 5.0f, 0.75f) },

        { "Pillar1",  4, 512, 512, "BakedPillar1.exr",
            BakeTargetType::Pillar, float3(-7.5f, 5.0f,  5.0f), float3(0.75f, 5.0f, 0.75f) },

        { "Pillar2",  5, 512, 512, "BakedPillar2.exr",
            BakeTargetType::Pillar, float3(7.5f, 5.0f, 15.0f), float3(0.75f, 5.0f, 0.75f) },

        { "Pillar3",  6, 512, 512, "BakedPillar3.exr",
            BakeTargetType::Pillar, float3(-7.5f, 5.0f, 15.0f), float3(0.75f, 5.0f, 0.75f) },

        { "Pillar4",  7, 512, 512, "BakedPillar4.exr",
            BakeTargetType::Pillar, float3(7.5f, 5.0f, 25.0f), float3(0.75f, 5.0f, 0.75f) },

        { "Pillar5",  8, 512, 512, "BakedPillar5.exr",
            BakeTargetType::Pillar, float3(-7.5f, 5.0f, 25.0f), float3(0.75f, 5.0f, 0.75f) },

        { "Pillar6",  9, 512, 512, "BakedPillar6.exr",
            BakeTargetType::Pillar, float3(7.5f, 5.0f, 35.0f), float3(0.75f, 5.0f, 0.75f) },

        { "Pillar7", 10, 512, 512, "BakedPillar7.exr",
            BakeTargetType::Pillar, float3(-7.5f, 5.0f, 35.0f), float3(0.75f, 5.0f, 0.75f) }
    };
}

void BakeLightMapSingle::loadDataSceneLightmaps()
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

void BakeLightMapSingle::setupCornellBakeTargets()
{
    mBakeTargets =
    {
        // Room shell targets (Quads)
        { "Floor",     0, 1024, 1024, "BakedCornell_Floor.exr"     },
        { "Ceiling",   1, 1024, 1024, "BakedCornell_Ceiling.exr"   },
        { "BackWall",  2, 1024, 1024, "BakedCornell_BackWall.exr"  },
        { "LeftWall",  3, 1024, 1024, "BakedCornell_LeftWall.exr"  },
        { "RightWall", 4, 1024, 1024, "BakedCornell_RightWall.exr" },

        // Light is instance ID 5; we skip baking the emissive light source.

        // Visibility-discontinuity test slab (Cube/Pillar)
        { "VisibilitySlab", 6, 512, 512, "BakedCornell_VisibilitySlab.exr",
            BakeTargetType::Pillar,
            float3(0.12f, 0.28f, 0.00f),       // translation
            float3(0.15f, 0.0175f, 0.11f),     // half extent (scaling * 0.5)
            float3(0.0f, 0.0f, 0.0f)           // rotation
        }
    };
}

void BakeLightMapSingle::loadCornellLightmaps()
{
    auto loadOne = [&](size_t idx, ref<Texture>& dst, const std::string& debugName)
        {
            if (idx >= mBakeTargets.size()) return;
            const auto& target = mBakeTargets[idx];
            dst = Texture::createFromFile(mpDevice, target.outputPath, true, false, ResourceBindFlags::ShaderResource);
            if (dst) dst->setName(debugName);
        };

    loadOne(0, mpCornellFloorLightmapShadowBoundaryTestScene, "CornellFloorLightmap");
    loadOne(1, mpCornellCeilingLightmapShadowBoundaryTestScene, "CornellCeilingLightmap");
    loadOne(2, mpCornellBackWallLightmapShadowBoundaryTestScene, "CornellBackWallLightmap");
    loadOne(3, mpCornellLeftWallLightmapShadowBoundaryTestScene, "CornellLeftWallLightmap");
    loadOne(4, mpCornellRightWallLightmapShadowBoundaryTestScene, "CornellRightWallLightmap");
    loadOne(5, mpCornellThinSlabLightmapShadowBoundaryTestScene, "CornellVisibilitySlabLightmap");
}

void BakeLightMapSingle::bindCornellData(ShaderVar applyVar)
{
    applyVar["gLinearSampler"] = mpLinearSampler;

    // Bind Textures
    applyVar["gFloorLightmap"] = mpCornellFloorLightmapShadowBoundaryTestScene;
    applyVar["gCeilingLightmap"] = mpCornellCeilingLightmapShadowBoundaryTestScene;
    applyVar["gBackWallLightmap"] = mpCornellBackWallLightmapShadowBoundaryTestScene;
    applyVar["gLeftWallLightmap"] = mpCornellLeftWallLightmapShadowBoundaryTestScene;
    applyVar["gRightWallLightmap"] = mpCornellRightWallLightmapShadowBoundaryTestScene;
    applyVar["gVisibilitySlabLightmap"] = mpCornellThinSlabLightmapShadowBoundaryTestScene;

    // Bind Instance IDs
    applyVar["PerFrameCB"]["gFloorInstanceID"] = 0;
    applyVar["PerFrameCB"]["gCeilingInstanceID"] = 1;
    applyVar["PerFrameCB"]["gBackWallInstanceID"] = 2;
    applyVar["PerFrameCB"]["gLeftWallInstanceID"] = 3;
    applyVar["PerFrameCB"]["gRightWallInstanceID"] = 4;
    applyVar["PerFrameCB"]["gVisibilitySlabInstanceID"] = 6;

    // Bind Pillar Data for Slab
    applyVar["PerFrameCB"]["gVisibilitySlabCenterW"] = float3(0.12f, 0.28f, 0.00f);
    applyVar["PerFrameCB"]["gVisibilitySlabHalfExtentW"] = float3(0.15f, 0.0175f, 0.11f);
}
