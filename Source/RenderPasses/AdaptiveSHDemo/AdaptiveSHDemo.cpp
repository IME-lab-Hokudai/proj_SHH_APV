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
#include "AdaptiveSHDemo.h"
#include "Rendering/Lights/EmissivePowerSampler.h"
#include "Rendering/Lights/EmissiveUniformSampler.h"
#include "Utils/Logger.h"

const char kStaticPassFile[] = "RenderPasses/AdaptiveSHDemo/StaticPass.slang";
const char kDynamicPassFile[] = "RenderPasses/AdaptiveSHDemo/DynamicPass.slang";

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
    return reflector;
}

void AdaptiveSHDemo::execute(RenderContext* pRenderContext, const RenderData& renderData)
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

        // ------------------------------------------------------------------
        // PASS 1: STATIC GEOMETRY (lightmaps)
        // ------------------------------------------------------------------
        auto applyVar = mpStaticVars->getRootVar();
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
        mpScene->rasterize(pRenderContext, mpGraphicsState.get(), mpStaticVars.get(), mpRasterState, mpRasterState);

        // ------------------------------------------------------------------
// PASS 2: DYNAMIC OBJECTS
// ------------------------------------------------------------------
        mpGraphicsState->setProgram(mpDynamicProgram);

        auto dynVar = mpDynamicVars->getRootVar();
        dynVar["gLinearSampler"] = mpLinearSampler;
        dynVar["PerFrameCB"]["gFirstDynamicInstanceID"] = mFirstDynamicInstanceID;
        dynVar["PerFrameCB"]["sampleIndex"] = mCurrentSample;
        mCurrentSample++;

        if (mpEmissiveSampler)
        {
            mpEmissiveSampler->bindShaderData(dynVar["PerFrameCB"]["emissiveSampler"]);
        }

        // Reuse the same FBO/depth. Do NOT clear again.
        mpScene->bindShaderDataForRaytracing(pRenderContext, dynVar["gScene"]);
        mpScene->rasterize(pRenderContext, mpGraphicsState.get(), mpDynamicVars.get(), mpRasterState, mpRasterState);

        // Restore static program if you want the state to remain consistent.
        mpGraphicsState->setProgram(mpStaticProgram);
    }
}

void AdaptiveSHDemo::renderUI(Gui::Widgets& widget) {}

void AdaptiveSHDemo::loadLightmaps()
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

void AdaptiveSHDemo::setScene(RenderContext* pRenderContext, const ref<Scene>& pScene)
{
    mpScene = pScene;
    if (mpScene)
    {
        mBakeTargets =
        {
            { "Floor",     0, 1024, 1024, "BakedFloor.exr"     },
            { "LeftWall",  1, 1024,  512, "BakedLeftWall.exr"  },
            { "RightWall", 2, 1024,  512, "BakedRightWall.exr" },
            { "RoofLeft", 11, 1024,  512, "BakedRoofLeft.exr"  },
            { "RoofRight",12, 1024,  512, "BakedRoofRight.exr" },
            { "Pillar0",  3, 512, 512, "BakedPillar0.exr", BakeTargetType::Pillar, float3(7.5f, 5.0f,  5.0f), float3(0.75f, 5.0f, 0.75f) },
            { "Pillar1",  4, 512, 512, "BakedPillar1.exr", BakeTargetType::Pillar, float3(-7.5f, 5.0f,  5.0f), float3(0.75f, 5.0f, 0.75f) },
            { "Pillar2",  5, 512, 512, "BakedPillar2.exr", BakeTargetType::Pillar, float3(7.5f, 5.0f, 15.0f), float3(0.75f, 5.0f, 0.75f) },
            { "Pillar3",  6, 512, 512, "BakedPillar3.exr", BakeTargetType::Pillar, float3(-7.5f, 5.0f, 15.0f), float3(0.75f, 5.0f, 0.75f) },
            { "Pillar4",  7, 512, 512, "BakedPillar4.exr", BakeTargetType::Pillar, float3(7.5f, 5.0f, 25.0f), float3(0.75f, 5.0f, 0.75f) },
            { "Pillar5",  8, 512, 512, "BakedPillar5.exr", BakeTargetType::Pillar, float3(-7.5f, 5.0f, 25.0f), float3(0.75f, 5.0f, 0.75f) },
            { "Pillar6",  9, 512, 512, "BakedPillar6.exr", BakeTargetType::Pillar, float3(7.5f, 5.0f, 35.0f), float3(0.75f, 5.0f, 0.75f) },
            { "Pillar7", 10, 512, 512, "BakedPillar7.exr", BakeTargetType::Pillar, float3(-7.5f, 5.0f, 35.0f), float3(0.75f, 5.0f, 0.75f) }
        };
        ProgramDesc previewDesc;
        previewDesc.addShaderModules(mpScene->getShaderModules());
        previewDesc.addShaderLibrary(kStaticPassFile)
            .vsEntry("vsMain")
            .psEntry("psMain");

        mpStaticProgram = Program::create(mpDevice, previewDesc, mpScene->getSceneDefines());
        mpStaticVars = ProgramVars::create(mpDevice, mpStaticProgram->getReflector());

        // rasterizer state
        RasterizerState::Desc rasterDesc;
        rasterDesc.setFillMode(RasterizerState::FillMode::Solid);
        rasterDesc.setCullMode(RasterizerState::CullMode::Back);
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
       
        loadLightmaps();
  
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
