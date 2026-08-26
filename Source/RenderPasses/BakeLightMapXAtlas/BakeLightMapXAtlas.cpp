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

#include "BakeLightMapXAtlas.h"

#include "Rendering/Lights/EmissivePowerSampler.h"
#include "Rendering/Lights/EmissiveUniformSampler.h"
#include "Utils/Logger.h"
#include "xatlas.h"

#include <algorithm>
#include <map>
#include <unordered_set>

namespace
{
    const char kUVRasterFile[] =
        "RenderPasses/BakeLightMapXAtlas/XAtlasUVRaster.slang";

    const char kBakingFile[] =
        "RenderPasses/BakeLightMapXAtlas/LightmapBakingSingle.rt.slang";

    const char kExtractFile[] =
        "RenderPasses/BakeLightMapXAtlas/ExtractTexelsSingle.cs.slang";

    const char kNormalizeFile[] =
        "RenderPasses/BakeLightMapXAtlas/NormalizeLightmapSingle.cs.slang";
}

extern "C" FALCOR_API_EXPORT void registerPlugin(
    Falcor::PluginRegistry& registry
)
{
    registry.registerClass<RenderPass, BakeLightMapXAtlas>();
}

BakeLightMapXAtlas::BakeLightMapXAtlas(
    ref<Device> pDevice,
    const Properties& props
)
    : RenderPass(pDevice)
{}

Properties BakeLightMapXAtlas::getProperties() const
{
    return {};
}

RenderPassReflection BakeLightMapXAtlas::reflect(
    const CompileData& compileData
)
{
    RenderPassReflection reflector;

    const uint2 sz = RenderPassHelpers::calculateIOSize(
        mOutputSizeSelection,
        mFixedOutputSize,
        compileData.defaultTexDims
    );

    reflector.addOutput("output", "Color")
        .texture2D(sz.x, sz.y, 4)
        .format(ResourceFormat::RGBA32Float);

    reflector.addOutput("depth", "Depth buffer")
        .format(ResourceFormat::D32Float)
        .bindFlags(ResourceBindFlags::DepthStencil)
        .texture2D(sz.x, sz.y, 4);

    return reflector;
}

void BakeLightMapXAtlas::execute(
    RenderContext* pRenderContext,
    const RenderData& renderData
)
{
    // Baking is currently triggered once from setScene(). Keeping execute()
    // empty avoids re-running the expensive atlas build every frame.
}

void BakeLightMapXAtlas::renderUI(Gui::Widgets& widget)
{}

void BakeLightMapXAtlas::setScene(
    RenderContext* pRenderContext,
    const ref<Scene>& pScene
)
{
    resetBakingState();

    mpScene = pScene;
    if (!mpScene)
        return;

    createUVRasterProgram();
    createComputePasses();
    createRayTracingProgram(pRenderContext);

    const std::vector<uint32_t> allInstanceIDs =
        collectTriangleInstanceIDs();

    if (allInstanceIDs.empty())
    {
        logWarning("BakeLightMapXAtlas: scene has no triangle-mesh instances.");
        return;
    }

    // Keep the current 48-instance Bistro test while the multi-page path is
    // being validated. Set mTestInstanceCount >= allInstanceIDs.size() to
    // build the complete scene atlas.
    const uint32_t instanceCount =
        std::min<uint32_t>(
            mTestInstanceCount,
            uint32_t(allInstanceIDs.size())
        );

    std::vector<uint32_t> instanceIDs(
        allInstanceIDs.begin(),
        allInstanceIDs.begin() + instanceCount
    );

    buildAtlasPages(
        pRenderContext,
        instanceIDs,
        mAtlasResolution,
        mTexelsPerUnit
    );

    logInfo(
        "Built {} lightmap atlas page(s) from {} instance(s).",
        mAtlasPages.size(),
        instanceIDs.size()
    );

    // Bake each physical xatlas page independently. All pages share the same
    // dimensions, exactly like the official xatlas multi-atlas example.
    //for (auto& page : mAtlasPages)
    //{
    //    bakeAtlasPage(
    //        pRenderContext,
    //        page
    //    );
    //}
}

void BakeLightMapXAtlas::resetBakingState()
{
    mAtlasPages.clear();

    mLightmapWidth = 0;
    mLightmapHeight = 0;
    mNumExtractedTexels = 0;
    mCurrentSample = 0;

    mpUVFbo = nullptr;
    mpTexelBuffer = nullptr;
    mpCounterBuffer = nullptr;
    mpAccumBuffer = nullptr;
    mpResultTex = nullptr;

    mpUVProgram = nullptr;
    mpUVVars = nullptr;
    mpUVGraphicsState = nullptr;
    mpRasterState = nullptr;
    mpProceduralVao = nullptr;

    mpExtractPass = nullptr;
    mpNormalizePass = nullptr;

    mpRtVars = nullptr;
    mpRtProgram = nullptr;

    mpEmissiveSampler.reset();
}

void BakeLightMapXAtlas::createUVRasterProgram()
{
    ProgramDesc desc;
    desc.addShaderModules(mpScene->getShaderModules());
    desc.addShaderLibrary(kUVRasterFile)
        .vsEntry("vsMain")
        .psEntry("psMain");

    mpUVProgram = Program::create(
        mpDevice,
        desc,
        mpScene->getSceneDefines()
    );

    mpUVVars = ProgramVars::create(
        mpDevice,
        mpUVProgram->getReflector()
    );

    RasterizerState::Desc rasterDesc;
    rasterDesc.setFillMode(RasterizerState::FillMode::Solid);
    rasterDesc.setCullMode(RasterizerState::CullMode::None);
    mpRasterState = RasterizerState::create(rasterDesc);

    DepthStencilState::Desc depthDesc;
    depthDesc.setDepthEnabled(false);
    ref<DepthStencilState> pDepthState =
        DepthStencilState::create(depthDesc);

    mpUVGraphicsState = GraphicsState::create(mpDevice);
    mpUVGraphicsState->setProgram(mpUVProgram);
    mpUVGraphicsState->setRasterizerState(mpRasterState);
    mpUVGraphicsState->setDepthStencilState(pDepthState);

    // The VS fetches the source geometry through gScene and uses only
    // SV_VertexID, so no scene vertex/index buffers are bound to the IA.
    mpProceduralVao = Vao::create(Vao::Topology::TriangleList);
    mpUVGraphicsState->setVao(mpProceduralVao);
}

void BakeLightMapXAtlas::createComputePasses()
{
    mpExtractPass = ComputePass::create(
        mpDevice,
        kExtractFile,
        "main",
        mpScene->getSceneDefines()
    );

    mpNormalizePass = ComputePass::create(
        mpDevice,
        kNormalizeFile,
        "main",
        mpScene->getSceneDefines()
    );
}

void BakeLightMapXAtlas::createRayTracingProgram(
    RenderContext* pRenderContext
)
{
    ProgramDesc rtDesc;
    rtDesc.addShaderModules(mpScene->getShaderModules());
    rtDesc.addShaderLibrary(kBakingFile);
    rtDesc.setMaxTraceRecursionDepth(3);
    rtDesc.setMaxPayloadSize(128);
    rtDesc.setMaxAttributeSize(8);
    rtDesc.addTypeConformances(mpScene->getTypeConformances());

    ref<RtBindingTable> sbt = RtBindingTable::create(
        2,
        2,
        mpScene->getGeometryCount()
    );

    sbt->setRayGen(rtDesc.addRayGen("rayGen"));
    sbt->setMiss(0, rtDesc.addMiss("primaryMiss"));
    sbt->setMiss(1, rtDesc.addMiss("shadowMiss"));

    auto primaryHit =
        rtDesc.addHitGroup("primaryClosestHit");
    auto shadowHit =
        rtDesc.addHitGroup("", "shadowAnyHit");

    const auto triangleGeometryIDs =
        mpScene->getGeometryIDs(Scene::GeometryType::TriangleMesh);

    sbt->setHitGroup(
        0,
        triangleGeometryIDs,
        primaryHit
    );

    sbt->setHitGroup(
        1,
        triangleGeometryIDs,
        shadowHit
    );

    const auto& pLights =
        mpScene->getILightCollection(pRenderContext);

    if (mpScene->useEmissiveLights())
    {
        FALCOR_ASSERT(
            pLights &&
            pLights->getActiveLightCount(pRenderContext) > 0
        );

        switch (mEmissiveSamplerType)
        {
        case EmissiveLightSamplerType::Uniform:
            mpEmissiveSampler =
                std::make_unique<EmissiveUniformSampler>(
                    pRenderContext,
                    pLights
                );
            break;

        case EmissiveLightSamplerType::LightBVH:
            mpEmissiveSampler =
                std::make_unique<LightBVHSampler>(
                    pRenderContext,
                    pLights,
                    mLightBVHOptions
                );
            break;

        case EmissiveLightSamplerType::Power:
            mpEmissiveSampler =
                std::make_unique<EmissivePowerSampler>(
                    pRenderContext,
                    pLights
                );
            break;

        default:
            FALCOR_THROW("Unknown emissive light sampler type.");
        }
    }

    mpRtProgram = Program::create(
        mpDevice,
        rtDesc,
        mpScene->getSceneDefines()
    );

    if (mpEmissiveSampler)
    {
        mpRtProgram->addDefines(
            mpEmissiveSampler->getDefines()
        );
    }

    DefineList lightDefines;
    lightDefines.add(
        "USE_ANALYTIC_LIGHTS",
        mpScene->useAnalyticLights() ? "1" : "0"
    );
    lightDefines.add(
        "USE_EMISSIVE_LIGHTS",
        mpScene->useEmissiveLights() ? "1" : "0"
    );

    mpRtProgram->addDefines(lightDefines);

    mpRtVars = RtProgramVars::create(
        mpDevice,
        mpRtProgram,
        sbt
    );
}

std::vector<uint32_t>
BakeLightMapXAtlas::collectTriangleInstanceIDs() const
{
    const auto instanceIDs =
        mpScene->getGeometryInstanceIDsByType(
            Scene::GeometryType::TriangleMesh
        );

    std::unordered_set<uint32_t> uniqueMeshes;
    uint64_t totalTriangles = 0;

    for (uint32_t instanceID : instanceIDs)
    {
        const auto& instance =
            mpScene->getGeometryInstance(instanceID);

        const MeshID meshID(instance.geometryID);
        const auto& mesh = mpScene->getMesh(meshID);

        uniqueMeshes.insert(meshID.get());
        totalTriangles += mesh.indexCount / 3;
    }

    logInfo(
        "Bistro lightmap targets: "
        "instances={} uniqueMeshes={} totalInstanceTriangles={}",
        instanceIDs.size(),
        uniqueMeshes.size(),
        totalTriangles
    );

    return instanceIDs;
}

void BakeLightMapXAtlas::buildAtlasPages(
    RenderContext* pRenderContext,
    const std::vector<uint32_t>& instanceIDs,
    uint32_t resolution,
    float texelsPerUnit
)
{
    if (instanceIDs.empty())
        FALCOR_THROW("Cannot build xatlas pages from an empty instance list.");

    if (resolution == 0)
        FALCOR_THROW("xatlas page resolution must be greater than zero.");

    if (texelsPerUnit <= 0.f)
        FALCOR_THROW("xatlas texelsPerUnit must be greater than zero.");

    struct InputMeshInfo
    {
        uint32_t instanceID = 0;
        MeshID meshID;
        std::vector<uint3> triangles;
    };

    std::vector<InputMeshInfo> inputMeshes;
    inputMeshes.reserve(instanceIDs.size());

    xatlas::Atlas* atlas = xatlas::Create();
    if (!atlas)
        FALCOR_THROW("xatlas::Create() failed.");

    // -----------------------------------------------------------------
    // Add one xatlas mesh per Falcor geometry instance.
    //
    // This mirrors the official xatlas example: atlas->meshes[i]
    // corresponds to the i-th successful AddMesh() call.
    // -----------------------------------------------------------------
    for (uint32_t instanceID : instanceIDs)
    {
        const auto& instance =
            mpScene->getGeometryInstance(instanceID);

        const MeshID meshID(instance.geometryID);
        const auto& mesh = mpScene->getMesh(meshID);

        const uint32_t vertexCount = mesh.vertexCount;
        const uint32_t triangleCount = mesh.indexCount / 3;

        ref<Buffer> pPositions =
            mpDevice->createStructuredBuffer(
                sizeof(float3),
                vertexCount,
                ResourceBindFlags::ShaderResource |
                ResourceBindFlags::UnorderedAccess
            );

        // Falcor's mesh extraction path expects this output even though
        // xatlas does not need the material UVs for this parameterization.
        ref<Buffer> pTexCoords =
            mpDevice->createStructuredBuffer(
                sizeof(float3),
                vertexCount,
                ResourceBindFlags::ShaderResource |
                ResourceBindFlags::UnorderedAccess
            );

        ref<Buffer> pTriangleIndices =
            mpDevice->createStructuredBuffer(
                sizeof(uint3),
                triangleCount,
                ResourceBindFlags::ShaderResource |
                ResourceBindFlags::UnorderedAccess
            );

        std::map<std::string, ref<Buffer>> buffers;
        buffers["positions"] = pPositions;
        buffers["texcrds"] = pTexCoords;
        buffers["triangleIndices"] = pTriangleIndices;

        mpScene->getMeshVerticesAndIndices(
            meshID,
            buffers
        );

        std::vector<float3> positions(vertexCount);
        std::vector<uint3> triangles(triangleCount);

        pPositions->getBlob(
            positions.data(),
            0,
            positions.size() * sizeof(float3)
        );

        pTriangleIndices->getBlob(
            triangles.data(),
            0,
            triangles.size() * sizeof(uint3)
        );

        std::vector<uint32_t> indices;
        indices.reserve(triangleCount * 3);

        for (const uint3& tri : triangles)
        {
            indices.push_back(tri.x);
            indices.push_back(tri.y);
            indices.push_back(tri.z);
        }

        xatlas::MeshDecl meshDecl{};
        meshDecl.vertexCount = vertexCount;
        meshDecl.vertexPositionData = positions.data();
        meshDecl.vertexPositionStride = sizeof(float3);
        meshDecl.indexCount = uint32_t(indices.size());
        meshDecl.indexData = indices.data();
        meshDecl.indexFormat = xatlas::IndexFormat::UInt32;

        const xatlas::AddMeshError error =
            xatlas::AddMesh(
                atlas,
                meshDecl,
                uint32_t(instanceIDs.size())
            );

        if (error != xatlas::AddMeshError::Success)
        {
            const std::string errorString =
                xatlas::StringForEnum(error);

            xatlas::Destroy(atlas);

            FALCOR_THROW(
                "xatlas::AddMesh() failed for instanceID={} meshID={}: {}",
                instanceID,
                meshID,
                errorString
            );
        }

        InputMeshInfo input;
        input.instanceID = instanceID;
        input.meshID = meshID;
        input.triangles = std::move(triangles);
        inputMeshes.push_back(std::move(input));

        logInfo(
            "xatlas input: instanceID={} meshID={} vertices={} triangles={}",
            instanceID,
            meshID,
            vertexCount,
            triangleCount
        );
    }

    // The official example calls this only to join any asynchronous AddMesh
    // work before printing totals. It is harmless here and makes the mesh
    // list complete before generation/validation.
    xatlas::AddMeshJoin(atlas);

    xatlas::ChartOptions chartOptions{};
    xatlas::PackOptions packOptions{};

    packOptions.resolution = resolution;
    packOptions.texelsPerUnit = texelsPerUnit;
    packOptions.padding = 4;
    packOptions.bilinear = true;
    packOptions.createImage = false;

    xatlas::Generate(
        atlas,
        chartOptions,
        packOptions
    );

    // NOW meshCount is valid.
    if (atlas->meshCount != inputMeshes.size())
    {
        const uint32_t outputMeshCount = atlas->meshCount;
        const size_t expectedMeshCount = inputMeshes.size();

        xatlas::Destroy(atlas);

        FALCOR_THROW(
            "xatlas mesh count mismatch after Generate(). Expected {}, got {}.",
            expectedMeshCount,
            outputMeshCount
        );
    }

    logInfo(
        "xatlas result: meshes={} charts={} pages={} size={}x{} "
        "texelsPerUnit={}",
        atlas->meshCount,
        atlas->chartCount,
        atlas->atlasCount,
        atlas->width,
        atlas->height,
        atlas->texelsPerUnit
    );

    if (atlas->atlasCount == 0 ||
        atlas->width == 0 ||
        atlas->height == 0)
    {
        xatlas::Destroy(atlas);
        FALCOR_THROW("xatlas generated no valid atlas pages.");
    }

    for (uint32_t pageIndex = 0;
        pageIndex < atlas->atlasCount;
        ++pageIndex)
    {
        logInfo(
            "xatlas page {} utilization={:.2f}%",
            pageIndex,
            atlas->utilization[pageIndex] * 100.f
        );
    }

    // -----------------------------------------------------------------
    // One AtlasPageData per xatlas atlasIndex, exactly like the example's
    // outputTrisImage[atlasIndex * imageDataSize].
    // -----------------------------------------------------------------
    mAtlasPages.clear();
    mAtlasPages.resize(atlas->atlasCount);

    for (uint32_t pageIndex = 0;
        pageIndex < atlas->atlasCount;
        ++pageIndex)
    {
        AtlasPageData& page = mAtlasPages[pageIndex];
        page.pageIndex = pageIndex;
        page.width = atlas->width;
        page.height = atlas->height;
        page.outputPath =
            "Bistro_AtlasPage_" +
            std::to_string(pageIndex) +
            ".exr";
    }

    uint64_t packedTriangleCount = 0;
    uint64_t unatlasedTriangleCount = 0;

    // -----------------------------------------------------------------
    // Convert xatlas output face-by-face.
    //
    // This deliberately follows the official example's triangle raster
    // loop:
    //   face f -> indexArray[f*3 + corner] -> Vertex
    //   Vertex::atlasIndex selects the physical page
    //   Vertex::uv is in atlas texel coordinates
    //
    // The example states atlasIndex is the same for all vertices in a face
    // and skips faces whose atlasIndex < 0.
    // -----------------------------------------------------------------
    for (uint32_t meshIndex = 0;
        meshIndex < atlas->meshCount;
        ++meshIndex)
    {
        const xatlas::Mesh& xaMesh =
            atlas->meshes[meshIndex];

        const InputMeshInfo& input =
            inputMeshes[meshIndex];

        const uint32_t triangleCount =
            uint32_t(input.triangles.size());

        if (xaMesh.indexCount != triangleCount * 3)
        {
            const uint32_t outputIndexCount = xaMesh.indexCount;
            xatlas::Destroy(atlas);

            FALCOR_THROW(
                "xatlas index-count mismatch for instanceID={}. "
                "Expected {}, got {}.",
                input.instanceID,
                triangleCount * 3,
                outputIndexCount
            );
        }

        // Index into page.instances for this input instance on each page.
        // -1 means the instance has not contributed a packed face to that
        // page yet.
        std::vector<int32_t> pageInstanceIndices(
            atlas->atlasCount,
            -1
        );

        for (uint32_t triangleID = 0;
            triangleID < triangleCount;
            ++triangleID)
        {
            const uint32_t firstIndex = triangleID * 3;

            const uint32_t xaIndex0 =
                xaMesh.indexArray[firstIndex + 0];
            const uint32_t xaIndex1 =
                xaMesh.indexArray[firstIndex + 1];
            const uint32_t xaIndex2 =
                xaMesh.indexArray[firstIndex + 2];

            const xatlas::Vertex& v0 =
                xaMesh.vertexArray[xaIndex0];
            const xatlas::Vertex& v1 =
                xaMesh.vertexArray[xaIndex1];
            const xatlas::Vertex& v2 =
                xaMesh.vertexArray[xaIndex2];

            // Preserve the face/corner correspondence we already validated:
            // xref maps each seam-split xatlas output vertex back to the
            // original input vertex.
            const uint3& sourceTriangle =
                input.triangles[triangleID];

            if (v0.xref != sourceTriangle.x ||
                v1.xref != sourceTriangle.y ||
                v2.xref != sourceTriangle.z)
            {
                xatlas::Destroy(atlas);

                FALCOR_THROW(
                    "xatlas face mapping mismatch for instanceID={} "
                    "triangleID={}: Falcor=({}, {}, {}) "
                    "xatlas xref=({}, {}, {}).",
                    input.instanceID,
                    triangleID,
                    sourceTriangle.x,
                    sourceTriangle.y,
                    sourceTriangle.z,
                    v0.xref,
                    v1.xref,
                    v2.xref
                );
            }

            const int32_t atlasIndex = v0.atlasIndex;

            // Official xatlas example: atlasIndex is identical for every
            // vertex in a face. Treat a violation as corrupted output.
            if (v1.atlasIndex != atlasIndex ||
                v2.atlasIndex != atlasIndex)
            {
                xatlas::Destroy(atlas);

                FALCOR_THROW(
                    "xatlas produced a cross-page face for instanceID={} "
                    "triangleID={}: ({}, {}, {}).",
                    input.instanceID,
                    triangleID,
                    v0.atlasIndex,
                    v1.atlasIndex,
                    v2.atlasIndex
                );
            }

            // Same behavior as example.cpp: skip faces that xatlas did not
            // atlas (for Bistro this includes the nearly-degenerate face we
            // already diagnosed).
            if (atlasIndex < 0)
            {
                ++unatlasedTriangleCount;
                continue;
            }

            if (atlasIndex >= int32_t(atlas->atlasCount))
            {
                xatlas::Destroy(atlas);

                FALCOR_THROW(
                    "xatlas returned invalid atlasIndex={} for "
                    "instanceID={} triangleID={}.",
                    atlasIndex,
                    input.instanceID,
                    triangleID
                );
            }

            AtlasPageData& page =
                mAtlasPages[uint32_t(atlasIndex)];

            int32_t& pageInstanceIndex =
                pageInstanceIndices[uint32_t(atlasIndex)];

            if (pageInstanceIndex < 0)
            {
                AtlasPageInstanceData instanceData;
                instanceData.instanceID = input.instanceID;
                instanceData.meshID = input.meshID;

                page.instances.push_back(
                    std::move(instanceData)
                );

                pageInstanceIndex =
                    int32_t(page.instances.size() - 1);
            }

            AtlasPageInstanceData& instanceData =
                page.instances[size_t(pageInstanceIndex)];

            TriangleLightmapUV uv;
            uv.uv0 = float2(
                v0.uv[0] / float(atlas->width),
                v0.uv[1] / float(atlas->height)
            );
            uv.uv1 = float2(
                v1.uv[0] / float(atlas->width),
                v1.uv[1] / float(atlas->height)
            );
            uv.uv2 = float2(
                v2.uv[0] / float(atlas->width),
                v2.uv[1] / float(atlas->height)
            );

            instanceData.triangleIDs.push_back(triangleID);
            instanceData.triangleUVs.push_back(uv);

            ++packedTriangleCount;
        }
    }

    xatlas::Destroy(atlas);

    createAtlasPageGpuBuffers();

    uint64_t storedTriangleCount = 0;

    for (const AtlasPageData& page : mAtlasPages)
    {
        uint64_t pageTriangleCount = 0;

        for (const AtlasPageInstanceData& instance : page.instances)
        {
            pageTriangleCount += instance.triangleIDs.size();
        }

        storedTriangleCount += pageTriangleCount;

        logInfo(
            "Atlas page {} stored: size={}x{} instances={} triangles={} "
            "output='{}'.",
            page.pageIndex,
            page.width,
            page.height,
            page.instances.size(),
            pageTriangleCount,
            page.outputPath
        );
    }

    if (storedTriangleCount != packedTriangleCount)
    {
        FALCOR_THROW(
            "Atlas storage mismatch. Packed {} triangles but stored {}.",
            packedTriangleCount,
            storedTriangleCount
        );
    }

    logInfo(
        "Atlas conversion complete: packedTriangles={} "
        "unatlasedTriangles={} pages={}.",
        packedTriangleCount,
        unatlasedTriangleCount,
        mAtlasPages.size()
    );
}

void BakeLightMapXAtlas::createAtlasPageGpuBuffers()
{
    for (AtlasPageData& page : mAtlasPages)
    {
        for (AtlasPageInstanceData& instance : page.instances)
        {
            if (instance.triangleIDs.size() !=
                instance.triangleUVs.size())
            {
                FALCOR_THROW(
                    "Atlas page {} instance {} has mismatched triangle "
                    "ID/UV counts ({} vs {}).",
                    page.pageIndex,
                    instance.instanceID,
                    instance.triangleIDs.size(),
                    instance.triangleUVs.size()
                );
            }

            const uint32_t triangleCount =
                uint32_t(instance.triangleIDs.size());

            if (triangleCount == 0)
                continue;

            instance.pTriangleIDBuffer =
                mpDevice->createStructuredBuffer(
                    sizeof(uint32_t),
                    triangleCount,
                    ResourceBindFlags::ShaderResource
                );

            instance.pTriangleIDBuffer->setBlob(
                instance.triangleIDs.data(),
                0,
                instance.triangleIDs.size() * sizeof(uint32_t)
            );

            instance.pTriangleUVBuffer =
                mpDevice->createStructuredBuffer(
                    sizeof(TriangleLightmapUV),
                    triangleCount,
                    ResourceBindFlags::ShaderResource
                );

            instance.pTriangleUVBuffer->setBlob(
                instance.triangleUVs.data(),
                0,
                instance.triangleUVs.size() * sizeof(TriangleLightmapUV)
            );
        }
    }
}

void BakeLightMapXAtlas::bakeAtlasPage(
    RenderContext* pRenderContext,
    AtlasPageData& page
)
{
    if (page.width == 0 || page.height == 0)
    {
        FALCOR_THROW(
            "Atlas page {} has invalid dimensions {}x{}.",
            page.pageIndex,
            page.width,
            page.height
        );
    }

    mLightmapWidth = page.width;
    mLightmapHeight = page.height;

    // -----------------------------------------------------------------
    // UV-space G-buffer for this one physical xatlas page.
    // -----------------------------------------------------------------
    ref<Texture> pPosTex =
        mpDevice->createTexture2D(
            mLightmapWidth,
            mLightmapHeight,
            ResourceFormat::RGBA32Float,
            1,
            1,
            nullptr,
            ResourceBindFlags::RenderTarget |
            ResourceBindFlags::ShaderResource
        );

    ref<Texture> pNormTex =
        mpDevice->createTexture2D(
            mLightmapWidth,
            mLightmapHeight,
            ResourceFormat::RGBA32Float,
            1,
            1,
            nullptr,
            ResourceBindFlags::RenderTarget |
            ResourceBindFlags::ShaderResource
        );

    pPosTex->setName(
        "XAtlasPage_Pos_" +
        std::to_string(page.pageIndex)
    );
    pNormTex->setName(
        "XAtlasPage_Normal_" +
        std::to_string(page.pageIndex)
    );

    mpUVFbo = Fbo::create(mpDevice);
    mpUVFbo->attachColorTarget(pPosTex, 0);
    mpUVFbo->attachColorTarget(pNormTex, 1);

    GraphicsState::Viewport viewport(
        0.f,
        0.f,
        float(mLightmapWidth),
        float(mLightmapHeight),
        0.f,
        1.f
    );

    mpUVGraphicsState->setViewport(0, viewport);
    mpUVGraphicsState->setFbo(mpUVFbo);

    pRenderContext->clearFbo(
        mpUVFbo.get(),
        float4(0.f),
        1.f,
        0,
        FboAttachmentType::Color
    );

    auto var = mpUVVars->getRootVar();
    mpScene->bindShaderData(var["gScene"]);

    uint64_t rasterTriangleCount = 0;

    // -----------------------------------------------------------------
    // Draw only the compact face subset assigned to this atlasIndex.
    //
    // XAtlasUVRaster.slang must use gTriangleIDs[localTriangleID] to map
    // the compact page-local triangle back to the original Falcor mesh
    // triangle before calling gScene.getIndices().
    // -----------------------------------------------------------------
    for (const AtlasPageInstanceData& instance : page.instances)
    {
        const uint32_t triangleCount =
            uint32_t(instance.triangleIDs.size());

        if (triangleCount == 0)
            continue;

        if (!instance.pTriangleIDBuffer ||
            !instance.pTriangleUVBuffer)
        {
            FALCOR_THROW(
                "Atlas page {} instance {} is missing GPU triangle data.",
                page.pageIndex,
                instance.instanceID
            );
        }

        var["gTriangleIDs"] =
            instance.pTriangleIDBuffer;

        var["gTriangleUVs"] =
            instance.pTriangleUVBuffer;

        var["BakeCB"]["gReceiverInstanceID"] =
            instance.instanceID;

        var["BakeCB"]["gTriangleCount"] =
            triangleCount;

        const uint32_t vertexCount =
            triangleCount * 3;

        pRenderContext->draw(
            mpUVGraphicsState.get(),
            mpUVVars.get(),
            vertexCount,
            0
        );

        rasterTriangleCount += triangleCount;
    }

    logInfo(
        "Rasterized atlas page {}: size={}x{} instances={} triangles={}.",
        page.pageIndex,
        mLightmapWidth,
        mLightmapHeight,
        page.instances.size(),
        rasterTriangleCount
    );

    // -----------------------------------------------------------------
    // Extract valid atlas texels.
    // -----------------------------------------------------------------
    const uint32_t totalTexels =
        mLightmapWidth * mLightmapHeight;

    mpTexelBuffer =
        mpDevice->createStructuredBuffer(
            sizeof(TexelSample),
            totalTexels
        );

    mpCounterBuffer =
        mpDevice->createBuffer(
            sizeof(uint32_t),
            ResourceBindFlags::UnorderedAccess,
            MemoryType::DeviceLocal
        );

    pRenderContext->clearUAV(
        mpCounterBuffer->getUAV().get(),
        uint4(0)
    );

    auto extractVar =
        mpExtractPass->getRootVar();

    extractVar["gPosW"] =
        mpUVFbo->getColorTexture(0);
    extractVar["gNormW"] =
        mpUVFbo->getColorTexture(1);
    extractVar["gTexelSamples"] =
        mpTexelBuffer;
    extractVar["gCounter"] =
        mpCounterBuffer;
    extractVar["CB"]["gResolution"] =
        uint2(mLightmapWidth, mLightmapHeight);

    mpExtractPass->execute(
        pRenderContext,
        mLightmapWidth,
        mLightmapHeight
    );

    uint32_t extractedCount = 0;
    mpCounterBuffer->getBlob(
        &extractedCount,
        0,
        sizeof(uint32_t)
    );

    mNumExtractedTexels = extractedCount;

    logInfo(
        "Atlas page {} extracted {} / {} texels ({:.2f}% coverage).",
        page.pageIndex,
        mNumExtractedTexels,
        totalTexels,
        totalTexels > 0
        ? 100.f * float(mNumExtractedTexels) / float(totalTexels)
        : 0.f
    );

    if (mNumExtractedTexels == 0)
    {
        logWarning(
            "Atlas page {} has no valid texels; skipping RT bake.",
            page.pageIndex
        );
        return;
    }

    // TexelSample::texelIndex is the original atlas pixel index, so the
    // accumulation buffer must contain one element per atlas pixel.
    mpAccumBuffer =
        mpDevice->createStructuredBuffer(
            sizeof(float4),
            totalTexels
        );

    pRenderContext->clearUAV(
        mpAccumBuffer->getUAV().get(),
        float4(0.f)
    );

    mCurrentSample = 0;

    for (uint32_t sample = 0;
        sample < mBakeSampleCount;
        ++sample)
    {
        traceOneSample(pRenderContext);
    }

    // -----------------------------------------------------------------
    // Normalize and save this physical page.
    // -----------------------------------------------------------------
    mpResultTex =
        mpDevice->createTexture2D(
            mLightmapWidth,
            mLightmapHeight,
            ResourceFormat::RGBA32Float,
            1,
            1,
            nullptr,
            ResourceBindFlags::ShaderResource |
            ResourceBindFlags::UnorderedAccess
        );

    pRenderContext->clearUAV(
        mpResultTex->getUAV().get(),
        float4(0.f)
    );

    auto normVar =
        mpNormalizePass->getRootVar();

    normVar["gAccumBuffer"] =
        mpAccumBuffer;
    normVar["gOutput"] =
        mpResultTex;
    normVar["CB"]["gTotalSamples"] =
        mCurrentSample;
    normVar["CB"]["gWidth"] =
        mLightmapWidth;

    mpNormalizePass->execute(
        pRenderContext,
        mLightmapWidth,
        mLightmapHeight
    );

    if (page.outputPath.empty())
    {
        FALCOR_THROW(
            "Atlas page {} has no output path.",
            page.pageIndex
        );
    }

    mpResultTex->captureToFile(
        0,
        0,
        page.outputPath,
        Bitmap::FileFormat::ExrFile,
        Bitmap::ExportFlags::Uncompressed
    );

    logInfo(
        "Saved atlas page {} to '{}' using {} samples.",
        page.pageIndex,
        page.outputPath,
        mCurrentSample
    );
}

void BakeLightMapXAtlas::traceOneSample(
    RenderContext* pRenderContext
)
{
    if (!mpTexelBuffer)
        FALCOR_THROW("Texel buffer is null.");

    if (!mpAccumBuffer)
        FALCOR_THROW("Accumulation buffer is null.");

    if (mNumExtractedTexels == 0)
        FALCOR_THROW("No extracted texels available for ray tracing.");

    auto rtVar = mpRtVars->getRootVar();

    rtVar["gTexelSamples"] =
        mpTexelBuffer;
    rtVar["gIrradianceAccum"] =
        mpAccumBuffer;

    rtVar["PerFrameCB"]["sampleIndex"] =
        mCurrentSample;
    rtVar["PerFrameCB"]["numTexels"] =
        mNumExtractedTexels;
    rtVar["PerFrameCB"]["bias"] =
        0.01f;

    if (mpEmissiveSampler)
    {
        mpEmissiveSampler->bindShaderData(
            rtVar["PerFrameCB"]["emissiveSampler"]
        );
    }

    mpScene->raytrace(
        pRenderContext,
        mpRtProgram.get(),
        mpRtVars,
        uint3(mNumExtractedTexels, 1, 1)
    );

    ++mCurrentSample;

    logInfo(
        "Completed RT sample {} for {} texels.",
        mCurrentSample,
        mNumExtractedTexels
    );
}
