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
#include "BakeLightMap.h"
#include "Rendering/Lights/EmissivePowerSampler.h"
#include "Rendering/Lights/EmissiveUniformSampler.h"
#include "BakeDataStructures.slang"

const char kBakingFile[] = "RenderPasses/BakeLightMap/LightmapBaking.rt.slang";
const char kUVrasterFile[] = "RenderPasses/BakeLightMap/UVPass.slang";
const char kExtractFile[] = "RenderPasses/BakeLightMap/ExtractTexels.cs.slang";
const char kNormalizeFile[] = "RenderPasses/BakeLightMap/NormalizeLightmap.cs.slang";
const char kShaderFile[] = "RenderPasses/BakeLightMap/ApplyLightmap.slang";

extern "C" FALCOR_API_EXPORT void registerPlugin(Falcor::PluginRegistry& registry)
{
    registry.registerClass<RenderPass, BakeLightMap>();
}

BakeLightMap::BakeLightMap(ref<Device> pDevice, const Properties& props) : RenderPass(pDevice)
{
    mpFbo = Fbo::create(mpDevice);
    Sampler::Desc samplerDesc;
    samplerDesc.setFilterMode(TextureFilteringMode::Linear, TextureFilteringMode::Linear, TextureFilteringMode::Linear);
    mpLinearSampler = mpDevice->createSampler(samplerDesc);
}

Properties BakeLightMap::getProperties() const
{
    return {};
}

RenderPassReflection BakeLightMap::reflect(const CompileData& compileData)
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

void BakeLightMap::execute(RenderContext* pRenderContext, const RenderData& renderData)
{
    // STEP 1: Rasterize Scene into UV Space
    auto pTargetFbo = renderData.getTexture("output");
    const float4 clearColor(0, 0, 0, 1);
    mpFbo->attachColorTarget(pTargetFbo, 0);
    // IMPORTANT: Set Viewport to match the 'output' texture size

    // Update frame dimension based on render pass output.
    auto pDepth = renderData.getTexture("depth");
    pRenderContext->clearDsv(pDepth->getDSV().get(), 1.f, 0);
    mpFbo->attachDepthStencilTarget(pDepth);
    pRenderContext->clearFbo(mpFbo.get(), clearColor, 1.0f, 0, FboAttachmentType::Color);
    if (mpScene) {
        if (mbloadLightMap) {


            auto applyVar = mpVars->getRootVar();
            applyVar["gLightmap"] = mpLoadedLightmap;
            applyVar["gLinearSampler"] = mpLinearSampler; // Standard linear sampler
            mpScene->rasterize(pRenderContext, mpGraphicsState.get(), mpVars.get(), mpRasterState, mpRasterState);
        }
        else {
            if (mNeedsPreparation) {
                pRenderContext->clearFbo(mpUVFbo.get(), float4(0), 1.0f, 0, FboAttachmentType::Color);

                mpScene->rasterize(pRenderContext, mpGraphicsState.get(), mpVars.get(), mpRasterState, mpRasterState);

                pRenderContext->clearUAV(mpCounterBuffer->getUAV().get(), uint4(0));
                auto var = mpExtractPass->getRootVar();
                var["gPosW"] = mpUVFbo->getColorTexture(0);
                var["gNormW"] = mpUVFbo->getColorTexture(1);
                var["gTexelSamples"] = mpTexelBuffer;
                var["gCounter"] = mpCounterBuffer;
                var["CB"]["gResolution"] = uint2(mLightmapWidth, mLightmapHeight);
                mpExtractPass->execute(pRenderContext, mLightmapWidth, mLightmapHeight);
                mNeedsPreparation = false; // Only do this again if scene changes
                mCurrentSample = 0;        // Reset baking progress

                // --- Readback the Texel Count ---
                uint32_t extractedCount = 0;
                // Since it is a 1-element buffer, we read 4 bytes (sizeof(uint32_t))
                mpCounterBuffer->getBlob(&extractedCount, 0, sizeof(uint32_t));

                // Store it in a class member so the RT pass can use it
                mNumExtractedTexels = extractedCount;
            }
            // -------------------------------------------------------------------------
            // STEP 3: ITERATIVE BAKING (Ray Tracing Accumulation)
            // ------------------------------------------------------------------------
            if (mCurrentSample < mTotalSamples)
            {
                if (mCurrentSample == 0) {
                    pRenderContext->clearUAV(mpAccumBuffer->getUAV().get(), float4(0));
                }
                auto rtVar = mpRtVars->getRootVar();

                // Bind the data we extracted in the preparation step
                rtVar["gTexelSamples"] = mpTexelBuffer;

                // Output and Sample Progress
                rtVar["gIrradianceAccum"] = mpAccumBuffer;
                rtVar["PerFrameCB"]["sampleIndex"] = mCurrentSample;
                rtVar["PerFrameCB"]["numTexels"] = mNumExtractedTexels;
                rtVar["PerFrameCB"]["bias"] = 0.01f;

                if (mpEmissiveSampler)
                    mpEmissiveSampler->bindShaderData(rtVar["PerFrameCB"]["emissiveSampler"]);

                // Launch Ray Tracer
                // We launch a 1D grid. The shader uses gCounter[0] to stop extra threads.
                uint32_t threadCount = mLightmapWidth * mLightmapHeight;
                mpScene->raytrace(pRenderContext, mpRtProgram.get(), mpRtVars, uint3(threadCount, 1, 1));

                mCurrentSample++;
                // STEP 4: FINALIZATION (Run only when we hit the target)
                if (mCurrentSample == mTotalSamples)
                {
                    auto var = mpNormalizePass->getRootVar();
                    // Bind the 1D Buffer and 2D Texture
                    var["gAccumBuffer"] = mpAccumBuffer;
                    var["gOutput"] = mpResultTex;

                    // Bind the constants
                    var["CB"]["gTotalSamples"] = mTotalSamples;
                    var["CB"]["gWidth"] = mLightmapWidth; // Essential for 2D -> 1D mapping

                    mpNormalizePass->execute(pRenderContext, mLightmapWidth, mLightmapHeight);

                    // --- NEW: Save to File ---
                    //ref<Texture> pOutputTex = renderData.getTexture("output");
                    if (mpResultTex)
                    {
                        // Define your path. You can use Falcor's Project Directory or a hardcoded one.
                        std::string path = "BakedLightmap.exr"; // Use .exr for high dynamic range
                        mpResultTex->captureToFile(0, 0, path, Bitmap::FileFormat::ExrFile, Bitmap::ExportFlags::Uncompressed);
                        printf("Lightmap saved to: %s\n", path.c_str());
                    }

                    printf("Lightmap baking complete (%u samples).\n", mTotalSamples);
                }
            }
        }
    }
}
void BakeLightMap::renderUI(Gui::Widgets& widget) {}

void BakeLightMap::loadLightmap(const std::filesystem::path& path)
{
    // Load as a 2D texture. Falcor handles EXR (HDR) automatically.
    // We set loadAsSrgb to false because lightmaps contain linear radiance data.
    //mpLoadedLightmap = Texture::createFromFile(mpDevice, path, true, false, ResourceBindFlags::ShaderResource, Bitmap::ImportFlags::ConvertToFloat16);
    mpLoadedLightmap = Texture::createFromFile(mpDevice, path, true, false, ResourceBindFlags::ShaderResource);
    mpLoadedLightmap->setName("lightMap");
    if (mpLoadedLightmap) {
        printf("Successfully loaded lightmap: %s\n", path.string().c_str());
    }
}

void BakeLightMap::setScene(RenderContext* pRenderContext, const ref<Scene>& pScene)
{
    mpScene = pScene;
    if (mpScene)
    {
        if (mbloadLightMap) {
            ProgramDesc previewDesc;
            previewDesc.addShaderModules(mpScene->getShaderModules());
            previewDesc.addShaderLibrary(kShaderFile)
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

            loadLightmap("BakedLightmap.exr");
        }
        else {
            // 1. Create the UV FBO (Position and Normal maps)
            ref<Texture> pPosTex = mpDevice->createTexture2D(mLightmapWidth, mLightmapHeight, ResourceFormat::RGBA32Float, 1, 1, nullptr, ResourceBindFlags::RenderTarget | ResourceBindFlags::ShaderResource);
            pPosTex->setName("PosTex");

            ref<Texture> pNormTex = mpDevice->createTexture2D(mLightmapWidth, mLightmapHeight, ResourceFormat::RGBA32Float, 1, 1, nullptr, ResourceBindFlags::RenderTarget | ResourceBindFlags::ShaderResource);
            pNormTex->setName("NormText");

            mpUVFbo = Fbo::create(mpDevice);
           
            mpUVFbo->attachColorTarget(pPosTex, 0);
            mpUVFbo->attachColorTarget(pNormTex, 1);
            uint32_t totalTexels = mLightmapWidth * mLightmapHeight;
            mpTexelBuffer = mpDevice->createStructuredBuffer(sizeof(TexelSample), totalTexels); // Adjust struct size to match TexelSample
            mpTexelBuffer->setName("TexelBuffer");
            mpCounterBuffer = mpDevice->createBuffer(sizeof(uint32_t), ResourceBindFlags::UnorderedAccess, MemoryType::DeviceLocal);
            mpCounterBuffer->setName("CounterBuffer");
            mpAccumBuffer = mpDevice->createStructuredBuffer(sizeof(float4), totalTexels);
            mpAccumBuffer->setName("AccumBuffer");
            mpResultTex = mpDevice->createTexture2D(mLightmapWidth, mLightmapHeight, ResourceFormat::RGBA16Float, 1, 1, nullptr, ResourceBindFlags::ShaderResource | ResourceBindFlags::RenderTarget | ResourceBindFlags::UnorderedAccess);
            mpResultTex->setName("BakeResultTex");
            ProgramDesc desc;
            desc.addShaderModules(mpScene->getShaderModules());
            desc.addShaderLibrary(kUVrasterFile)
                .vsEntry("vsMain")  // Vertex shader entry point
                .psEntry("psMain"); // Pixel shader entry point;
            mpProgram = Program::create(mpDevice, desc, mpScene->getSceneDefines());
            mpVars = ProgramVars::create(mpDevice, mpProgram->getReflector());

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
            mpGraphicsState->setProgram(mpProgram);
            mpGraphicsState->setRasterizerState(mpRasterState);
            mpGraphicsState->setFbo(mpUVFbo);
            mpGraphicsState->setDepthStencilState(pDsState);
            GraphicsState::Viewport vp(0.0f, 0.0f, (float)mLightmapWidth, (float)mLightmapHeight, 0.0f, 1.0f);
            mpGraphicsState->setViewport(0, vp);

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
