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
 ***************************************************************************/
#pragma once

#include "Falcor.h"
#include "RenderGraph/RenderPass.h"
#include "RenderGraph/RenderPassHelpers.h"
#include "Rendering/Lights/EmissiveLightSampler.h"
#include "Rendering/Lights/LightBVHSampler.h"
#include "BakeDataStructures.slang"

#include <filesystem>
#include <limits>
#include <unordered_map>
#include <vector>

using namespace Falcor;

class BakeLightMapXAtlas : public RenderPass
{
public:
    FALCOR_PLUGIN_CLASS(
        BakeLightMapXAtlas,
        "BakeLightMapXAtlas",
        "Bake scene lightmaps using xatlas-generated multi-page UV atlases."
    );

    static ref<BakeLightMapXAtlas> create(
        ref<Device> pDevice,
        const Properties& props
    )
    {
        return make_ref<BakeLightMapXAtlas>(pDevice, props);
    }

    BakeLightMapXAtlas(
        ref<Device> pDevice,
        const Properties& props
    );

    Properties getProperties() const override;
    RenderPassReflection reflect(const CompileData& compileData) override;
    void compile(RenderContext* pRenderContext, const CompileData& compileData) override {}
    void execute(RenderContext* pRenderContext, const RenderData& renderData) override;
    void renderUI(Gui::Widgets& widget) override;
    void setScene(RenderContext* pRenderContext, const ref<Scene>& pScene) override;
    bool onMouseEvent(const MouseEvent& mouseEvent) override { return false; }
    bool onKeyEvent(const KeyboardEvent& keyEvent) override { return false; }

private:
    struct AtlasPageTriangleData
    {
        uint32_t instanceID = 0;
        uint32_t meshID = 0;
        uint32_t triangleID = 0;
        TriangleLightmapUV uv;
    };

    struct AtlasPageData
    {
        uint32_t pageIndex = 0;
        uint32_t width = 0;
        uint32_t height = 0;
        std::filesystem::path outputPath;

        std::vector<AtlasPageTriangleData> triangles;

        ref<Buffer> pInstanceIDBuffer;
        ref<Buffer> pTriangleIDBuffer;
        ref<Buffer> pTriangleUVBuffer;
    };

    // CPU copy of one unique Falcor mesh. Reused instances share this geometry
    // during xatlas submission, avoiding repeated GPU extraction/readback.
    struct CachedMeshGeometry
    {
        std::vector<float3> positions;
        std::vector<uint3> triangles;
        std::vector<uint32_t> indices;
    };


    void resetBakingState();

    void createUVRasterProgram();
    void createComputePasses();
    void createRayTracingProgram(RenderContext* pRenderContext);

    std::vector<uint32_t> collectTriangleInstanceIDs() const;
    void buildMeshGeometryCache(const std::vector<uint32_t>& instanceIDs);

    // Build one xatlas object for the complete selected scene. This follows the
    // official xatlas examples: add all meshes, ComputeCharts once, then
    // PackCharts once. texelsPerUnit is left at 0 so xatlas estimates a density
    // that approximately matches targetResolution.
    std::vector<AtlasPageData> buildGlobalAtlas(
        const std::vector<uint32_t>& instanceIDs,
        uint32_t targetResolution,
        float& outChosenTexelsPerUnit
    );

    void buildAndSaveAtlasMapping(
        const std::vector<uint32_t>& instanceIDs
    );

    std::vector<AtlasPageData> loadAtlasPagesFromMapping() const;

    void saveAtlasManifest(
        uint32_t pageCount,
        uint32_t width,
        uint32_t height,
        uint64_t triangleRecordCount,
        float texelsPerUnit
    ) const;

    void createAtlasPageGpuBuffers(AtlasPageData& page);
    void releaseAtlasPageGpuBuffers(AtlasPageData& page);
    void bakeSavedAtlasPages(RenderContext* pRenderContext);

    void bakeAtlasPage(
        RenderContext* pRenderContext,
        AtlasPageData& page
    );

    void traceOneSample(RenderContext* pRenderContext);

private:
    ref<Scene> mpScene;

    ref<Program> mpUVProgram;
    ref<ProgramVars> mpUVVars;
    ref<GraphicsState> mpUVGraphicsState;
    ref<RasterizerState> mpRasterState;
    ref<Vao> mpProceduralVao;
    ref<Fbo> mpUVFbo;

    ref<ComputePass> mpExtractPass;
    ref<ComputePass> mpNormalizePass;

    ref<Program> mpRtProgram;
    ref<RtProgramVars> mpRtVars;

    EmissiveLightSamplerType mEmissiveSamplerType =
        EmissiveLightSamplerType::Uniform;
    std::unique_ptr<EmissiveLightSampler> mpEmissiveSampler;
    mutable LightBVHSampler::Options mLightBVHOptions;

    ref<Buffer> mpTexelBuffer;
    ref<Buffer> mpCounterBuffer;
    ref<Buffer> mpAccumBuffer;
    ref<Texture> mpResultTex;

    std::unordered_map<uint32_t, CachedMeshGeometry> mMeshGeometryCache;

    RenderPassHelpers::IOSize mOutputSizeSelection =
        RenderPassHelpers::IOSize::Default;
    uint2 mFixedOutputSize = { 512, 512 };

    uint32_t mLightmapWidth = 0;
    uint32_t mLightmapHeight = 0;
    uint32_t mNumExtractedTexels = 0;
    uint32_t mCurrentSample = 0;

    // Atlas/bake settings.
    //
    // xatlas packing mode:
    //   resolution      = mAtlasResolution
    //   texelsPerUnit   = 0
    //
    // Per xatlas documentation this asks xatlas to estimate the texel density
    // so the complete input approximately matches the requested resolution,
    // instead of enforcing a fixed density and creating many atlas pages.
    uint32_t mAtlasResolution = 2048;
    uint32_t mBakeSampleCount = 64;

    // Set false after a successful full atlas build to reuse
    // Bistro_AtlasMapping.bin and skip xatlas completely on later high-spp runs.
    bool mRebuildAtlas = true;

    // Temporary validation limit. Leave at max for the complete Bistro scene.
    uint32_t mTestInstanceCount = std::numeric_limits<uint32_t>::max();

    bool mBakeCompleted = false;
};
