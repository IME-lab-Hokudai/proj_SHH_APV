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
#pragma once
#include "Falcor.h"
#include "RenderGraph/RenderPass.h"
#include "RenderGraph/RenderPassHelpers.h"
#include "Rendering/Lights/EmissiveLightSampler.h"
#include "Rendering/Lights/LightBVHSampler.h"
#include "AdaptiveProbeVolume.h"
#include "UniformProbeVolume.h"
#include "ProbeVisualizePass.h"
using namespace Falcor;

class AdaptiveSHDemo : public RenderPass
{
public:
    FALCOR_PLUGIN_CLASS(AdaptiveSHDemo, "AdaptiveSHDemo", "Insert pass description here.");

    static ref<AdaptiveSHDemo> create(ref<Device> pDevice, const Properties& props)
    {
        return make_ref<AdaptiveSHDemo>(pDevice, props);
    }

    AdaptiveSHDemo(ref<Device> pDevice, const Properties& props);

    virtual Properties getProperties() const override;
    virtual RenderPassReflection reflect(const CompileData& compileData) override;
    virtual void compile(RenderContext* pRenderContext, const CompileData& compileData) override {}
    virtual void execute(RenderContext* pRenderContext, const RenderData& renderData) override;
    virtual void renderUI(Gui::Widgets& widget) override;
    void loadLightmaps();
    virtual void setScene(RenderContext* pRenderContext, const ref<Scene>& pScene) override;
    virtual bool onMouseEvent(const MouseEvent& mouseEvent) override { return false; }
    virtual bool onKeyEvent(const KeyboardEvent& keyEvent) override { return false; }
    void loadDataSceneLightmaps();
    void setupDataSceneBakeTargets();
    void bindDataSceneData(ShaderVar var);
private:
    enum class BakeTargetType
    {
        Quad,
        Pillar
    };

    struct BakeTarget
    {
        std::string name;
        uint32_t instanceID;
        uint32_t width;
        uint32_t height;
        std::string outputPath;
        BakeTargetType type = BakeTargetType::Quad;
        float3 pillarCenterW = float3(0.f);
        float3 pillarHalfExtentW = float3(1.f);
        float3 pillarRotationEulerDeg = float3(0.f);
    };

    ref<Scene> mpScene;
    ref<Program> mpStaticProgram;
    ref<ProgramVars> mpStaticVars;
    ref<GraphicsState> mpGraphicsState;
    ref<RasterizerState> mpRasterState;


    ref<Fbo> mpFbo;

    RenderPassHelpers::IOSize mOutputSizeSelection = RenderPassHelpers::IOSize::Default;
    uint2 mFixedOutputSize = { 512, 512 };
    EmissiveLightSamplerType mEmissiveSamplerType = EmissiveLightSamplerType::Uniform;
    std::unique_ptr<EmissiveLightSampler> mpEmissiveSampler;
    mutable LightBVHSampler::Options mLightBVHOptions;
    ref<Program> mpRtProgram;
    ref<RtProgramVars> mpRtVars;

    ref<Sampler> mpLinearSampler;

    std::vector<BakeTarget> mBakeTargets;
    uint32_t mCurrentTargetIndex = 0;

    //ref<Texture> mpFloorLightmap;
    //ref<Texture> mpLeftWallLightmap;
    //ref<Texture> mpRightWallLightmap;
    //ref<Texture> mpRoofLeftLightmap;
    //ref<Texture> mpRoofRightLightmap;
    //ref<Texture> mpPillar0Lightmap;
    //ref<Texture> mpPillar1Lightmap;
    //ref<Texture> mpPillar2Lightmap;
    //ref<Texture> mpPillar3Lightmap;
    //ref<Texture> mpPillar4Lightmap;
    //ref<Texture> mpPillar5Lightmap;
    //ref<Texture> mpPillar6Lightmap;
    //ref<Texture> mpPillar7Lightmap;

    //uint32_t mFloorInstanceID = 0;
    //uint32_t mLeftWallInstanceID = 1;
    //uint32_t mRightWallInstanceID = 2;
    //uint32_t mRoofLeftInstanceID = 11;
    //uint32_t mRoofRightInstanceID = 12;
    //uint32_t mPillar0InstanceID = 3;
    //uint32_t mPillar1InstanceID = 4;
    //uint32_t mPillar2InstanceID = 5;
    //uint32_t mPillar3InstanceID = 6;
    //uint32_t mPillar4InstanceID = 7;
    //uint32_t mPillar5InstanceID = 8;
    //uint32_t mPillar6InstanceID = 9;
    //uint32_t mPillar7InstanceID = 10;

    // First dynamic object for now.
    // Dynamic raster pass
    ref<Program> mpDynamicProgram;
    ref<ProgramVars> mpDynamicVars;
    uint32_t mCurrentSample = 0;
    uint32_t mFirstDynamicInstanceID = 27;
    ref<ComputePass> mpFilterPass;

    ref<Program> mpCompositeProgram;
    ref<ProgramVars> mpCompositeVars;
    ref<GraphicsState> mpCompositeState;
    ref<Vao> mpEmptyVao;

    ref<AdaptiveProbeVolume> mAdaptiveProbeVolume;
    ref<UniformProbeVolume> mUniformProbeVolume;

    //visualize probe grid
    ref<ProbeVisualizePass> mpProbeVisualizePass;
    bool mbShowAdaptiveGrid = false;
    // Add this array to track checkbox states
    bool mVisLevels[8] = { true, true, true, true, true, true, true, true };
    bool mbDrawLeafOnly = false;

    ref<Texture> mpDataFloorLightmap;
    ref<Texture> mpDataCeilingLightmap;
    ref<Texture> mpDataLeftWallLightmap;
    ref<Texture> mpDataRightWallLightmap;
    ref<Texture> mpDataBackWallLightmap;
    ref<Texture> mpDataFrontWallLightmap;

    ref<Texture> mpDataTallBoxALightmap;
    ref<Texture> mpDataWideBoxBLightmap;
    ref<Texture> mpDataBlockCLightmap;
    ref<Texture> mpDataLowBoxDLightmap;
    ref<Texture> mpDataThinSlabELightmap;
    ref<Texture> mpDataThinSlabFLightmap;
    ref<Texture> mpDataTallBoxJLightmap;
    ref<Texture> mpDataShortBoxLLightmap;
};
