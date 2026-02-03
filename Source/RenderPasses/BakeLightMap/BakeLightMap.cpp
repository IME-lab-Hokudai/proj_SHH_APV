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

const char kBakingFile[] = "RenderPasses/BakeLightMap/LightmapBaking.rt.slang";

extern "C" FALCOR_API_EXPORT void registerPlugin(Falcor::PluginRegistry& registry)
{
    registry.registerClass<RenderPass, BakeLightMap>();
}

BakeLightMap::BakeLightMap(ref<Device> pDevice, const Properties& props) : RenderPass(pDevice) {}

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
    return reflector;
}

void BakeLightMap::execute(RenderContext* pRenderContext, const RenderData& renderData)
{
    // renderData holds the requested resources
    // auto& pTexture = renderData.getTexture("src");
}

void BakeLightMap::renderUI(Gui::Widgets& widget) {}

void BakeLightMap::setScene(RenderContext* pRenderContext, const ref<Scene>& pScene)
{
    mpScene = pScene;
    if (mpScene)
    {
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
