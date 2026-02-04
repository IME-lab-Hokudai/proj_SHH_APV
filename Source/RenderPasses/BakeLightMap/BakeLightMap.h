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

using namespace Falcor;

class BakeLightMap : public RenderPass
{
public:
    FALCOR_PLUGIN_CLASS(BakeLightMap, "BakeLightMap", "Insert pass description here.");

    static ref<BakeLightMap> create(ref<Device> pDevice, const Properties& props)
    {
        return make_ref<BakeLightMap>(pDevice, props);
    }

    BakeLightMap(ref<Device> pDevice, const Properties& props);

    virtual Properties getProperties() const override;
    virtual RenderPassReflection reflect(const CompileData& compileData) override;
    virtual void compile(RenderContext* pRenderContext, const CompileData& compileData) override {}
    virtual void execute(RenderContext* pRenderContext, const RenderData& renderData) override;
    virtual void renderUI(Gui::Widgets& widget) override;
    void loadLightmap(const std::filesystem::path& path);
    virtual void setScene(RenderContext* pRenderContext, const ref<Scene>& pScene) override;
    virtual bool onMouseEvent(const MouseEvent& mouseEvent) override { return false; }
    virtual bool onKeyEvent(const KeyboardEvent& keyEvent) override { return false; }

private:

    ref<Scene> mpScene;
    ref<Program> mpProgram;
    ref<ProgramVars> mpVars;
    ref<GraphicsState> mpGraphicsState;
    ref<RasterizerState> mpRasterState;
    ref<Fbo> mpUVFbo;
    ref<Fbo> mpFbo;

    RenderPassHelpers::IOSize mOutputSizeSelection = RenderPassHelpers::IOSize::Default;
    uint2 mFixedOutputSize = { 512, 512 };
    uint32_t mLightmapWidth = 512;
    uint32_t mLightmapHeight = 512;
    EmissiveLightSamplerType mEmissiveSamplerType = EmissiveLightSamplerType::Uniform;
    std::unique_ptr<EmissiveLightSampler> mpEmissiveSampler;
    mutable LightBVHSampler::Options mLightBVHOptions;
    ref<Program> mpRtProgram;
    ref<RtProgramVars> mpRtVars;

    ref<ComputePass> mpExtractPass;
    ref<ComputePass> mpNormalizePass;
    ref<Buffer> mpTexelBuffer;
    ref<Buffer> mpCounterBuffer;
    ref<Buffer> mpAccumBuffer;
    // In BakeLightMap.h
    ref<Texture> mpResultTex;
    uint32_t mNumExtractedTexels = 0;
    bool mNeedsPreparation = true;
    uint32_t mCurrentSample = 0;
    uint32_t mTotalSamples = 64; // Set your target quality here
    bool mbloadLightMap = true;
    ref<Texture> mpLoadedLightmap;
    ref<Sampler> mpLinearSampler;
};
