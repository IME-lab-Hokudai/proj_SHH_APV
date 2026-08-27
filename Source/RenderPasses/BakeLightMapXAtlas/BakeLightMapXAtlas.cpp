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
#include <array>
#include <atomic>
#include <cstdint>
#include <cstring>
#include <fstream>
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

    const char kAtlasMappingFile[] = "Bistro_AtlasMapping.bin";
    const char kAtlasManifestFile[] = "Bistro_AtlasManifest.txt";

    //const char kAtlasMappingFile[] = "Room_AtlasMapping.bin";
    //const char kAtlasManifestFile[] = "Room_AtlasManifest.txt";
    struct XAtlasProgressState
    {
        std::array<std::atomic<int>, 4> lastReported;

        XAtlasProgressState()
        {
            for (auto& value : lastReported)
                value.store(-10);
        }
    };

    bool xatlasProgressCallback(
        xatlas::ProgressCategory category,
        int progress,
        void* pUserData
    )
    {
        auto* pState = static_cast<XAtlasProgressState*>(pUserData);
        const size_t categoryIndex = static_cast<size_t>(category);

        if (!pState || categoryIndex >= pState->lastReported.size())
            return true;

        int previous = pState->lastReported[categoryIndex].load();

        if (progress == previous)
            return true;

        // xatlas may call the callback from worker threads. Only log when the
        // category advances by at least 10%, plus the final 100% update.
        while (progress == 100 || progress >= previous + 10)
        {
            if (pState->lastReported[categoryIndex].compare_exchange_weak(
                previous,
                progress
            ))
            {
                const char* categoryName = "Unknown";

                switch (category)
                {
                case xatlas::ProgressCategory::AddMesh:
                    categoryName = "AddMesh";
                    break;
                case xatlas::ProgressCategory::ComputeCharts:
                    categoryName = "ComputeCharts";
                    break;
                case xatlas::ProgressCategory::PackCharts:
                    categoryName = "PackCharts";
                    break;
                case xatlas::ProgressCategory::BuildOutputMeshes:
                    categoryName = "BuildOutputMeshes";
                    break;
                }

                logInfo(
                    "xatlas progress: {} {}%",
                    categoryName,
                    progress
                );
                break;
            }
        }

        return true;
    }

#pragma pack(push, 1)
    struct AtlasMappingHeader
    {
        char magic[8] = { 'X', 'A', 'L', 'M', 'A', 'P', '0', '1' };
        uint32_t version = 1;
        uint32_t pageCount = 0;
        uint32_t width = 0;
        uint32_t height = 0;
        float texelsPerUnit = 0.f;
        uint64_t triangleRecordCount = 0;
    };

    struct AtlasMappingRecord
    {
        uint32_t pageIndex = 0;
        uint32_t instanceID = 0;
        uint32_t meshID = 0;
        uint32_t triangleID = 0;
        float uv[6] = {};
    };
#pragma pack(pop)

    static_assert(sizeof(AtlasMappingHeader) == 36);
    static_assert(sizeof(AtlasMappingRecord) == 40);
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
    if (!mpScene || mBakeCompleted)
        return;

    const bool mappingExists =
        std::filesystem::exists(kAtlasMappingFile);

    if (mRebuildAtlas || !mappingExists)
    {
        if (!mRebuildAtlas && !mappingExists)
        {
            logWarning(
                "Atlas reuse requested, but '{}' does not exist. "
                "Rebuilding the atlas mapping.",
                kAtlasMappingFile
            );
        }

        const std::vector<uint32_t> allInstanceIDs =
            collectTriangleInstanceIDs();

        if (allInstanceIDs.empty())
        {
            logWarning(
                "BakeLightMapXAtlas: scene has no triangle-mesh instances."
            );
            mBakeCompleted = true;
            return;
        }

        const uint32_t instanceCount =
            std::min<uint32_t>(
                mTestInstanceCount,
                uint32_t(allInstanceIDs.size())
            );

        std::vector<uint32_t> instanceIDs(
            allInstanceIDs.begin(),
            allInstanceIDs.begin() + instanceCount
        );

        logInfo(
            "Starting global xatlas build: instances={} / {} "
            "targetResolution={} automaticTexelDensity=true.",
            instanceIDs.size(),
            allInstanceIDs.size(),
            mAtlasResolution
        );

        // Build and finalize the complete persistent mapping BEFORE any RT
        // baking. If a later high-spp bake is interrupted, the UV work is
        // already safely reusable on the next run.
        buildAndSaveAtlasMapping(instanceIDs);
    }
    else
    {
        logInfo(
            "Reusing existing xatlas mapping '{}'. "
            "No xatlas charting/packing will run.",
            kAtlasMappingFile
        );
    }

    // GPU programs are only needed for the actual lightmap bake. Creating
    // them here also keeps setScene() lightweight and lets Falcor finish its
    // normal scene update before the one-shot work begins.
    if (!mpUVProgram)
        createUVRasterProgram();
    if (!mpExtractPass || !mpNormalizePass)
        createComputePasses();
    if (!mpRtProgram || !mpRtVars)
        createRayTracingProgram(pRenderContext);

    bakeSavedAtlasPages(pRenderContext);

    mBakeCompleted = true;
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
    mBakeCompleted = false;
    if (mpScene) {
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

void BakeLightMapXAtlas::resetBakingState()
{
    mMeshGeometryCache.clear();

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
        "instances={} uniqueMeshes={} totalInstanceTriangles={}.",
        instanceIDs.size(),
        uniqueMeshes.size(),
        totalTriangles
    );

    return instanceIDs;
}

void BakeLightMapXAtlas::buildMeshGeometryCache(
    const std::vector<uint32_t>& instanceIDs
)
{
    mMeshGeometryCache.clear();

    std::unordered_set<uint32_t> requiredMeshIDs;
    requiredMeshIDs.reserve(instanceIDs.size());

    for (uint32_t instanceID : instanceIDs)
    {
        const auto& instance =
            mpScene->getGeometryInstance(instanceID);
        requiredMeshIDs.insert(MeshID(instance.geometryID).get());
    }

    mMeshGeometryCache.reserve(requiredMeshIDs.size());

    uint32_t cacheIndex = 0;
    uint64_t cachedTriangles = 0;

    for (uint32_t meshIDValue : requiredMeshIDs)
    {
        const MeshID meshID(meshIDValue);
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

        // Falcor's extraction helper expects a texcoord output even though
        // xatlas only consumes positions/indices here.
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

        mpScene->getMeshVerticesAndIndices(meshID, buffers);

        CachedMeshGeometry geometry;
        geometry.positions.resize(vertexCount);
        geometry.triangles.resize(triangleCount);
        geometry.indices.reserve(size_t(triangleCount) * 3u);

        pPositions->getBlob(
            geometry.positions.data(),
            0,
            geometry.positions.size() * sizeof(float3)
        );

        pTriangleIndices->getBlob(
            geometry.triangles.data(),
            0,
            geometry.triangles.size() * sizeof(uint3)
        );

        for (const uint3& tri : geometry.triangles)
        {
            geometry.indices.push_back(tri.x);
            geometry.indices.push_back(tri.y);
            geometry.indices.push_back(tri.z);
        }

        cachedTriangles += triangleCount;
        mMeshGeometryCache.emplace(
            meshIDValue,
            std::move(geometry)
        );

        ++cacheIndex;
        if (cacheIndex == uint32_t(requiredMeshIDs.size()) ||
            (cacheIndex % 100u) == 0u)
        {
            logInfo(
                "Mesh cache: {}/{} unique meshes extracted.",
                cacheIndex,
                requiredMeshIDs.size()
            );
        }
    }

    logInfo(
        "Mesh cache complete: uniqueMeshes={} uniqueTriangles={}.",
        mMeshGeometryCache.size(),
        cachedTriangles
    );
}

std::vector<BakeLightMapXAtlas::AtlasPageData>
BakeLightMapXAtlas::buildGlobalAtlas(
    const std::vector<uint32_t>& instanceIDs,
    uint32_t targetResolution,
    float& outChosenTexelsPerUnit
)
{
    if (instanceIDs.empty())
        FALCOR_THROW("Cannot build xatlas from an empty instance list.");

    if (targetResolution == 0)
        FALCOR_THROW("xatlas target resolution must be greater than zero.");

    struct InputMeshInfo
    {
        uint32_t instanceID = 0;
        uint32_t meshID = 0;
        uint32_t triangleCount = 0;
    };

    std::vector<InputMeshInfo> inputMeshes;
    inputMeshes.reserve(instanceIDs.size());

    xatlas::Atlas* atlas = xatlas::Create();
    if (!atlas)
        FALCOR_THROW("xatlas::Create() failed.");

    XAtlasProgressState progressState;
    xatlas::SetProgressCallback(
        atlas,
        xatlasProgressCallback,
        &progressState
    );

    uint64_t submittedTriangles = 0;

    const AnimationController* pAnimationController =
        mpScene->getAnimationController();

    if (!pAnimationController)
    {
        xatlas::Destroy(atlas);
        FALCOR_THROW(
            "Scene has no AnimationController."
        );
    }

    const auto& globalMatrices =
        pAnimationController->getGlobalMatrices();

    // -----------------------------------------------------------------
    // Official-example-style submission:
    // one xatlas::Atlas containing every selected Falcor instance.
    //
    // We still submit one xatlas mesh per Falcor geometry instance because
    // different instances can receive different baked lighting even when they
    // reference the same source mesh. Cached CPU geometry avoids repeated
    // Falcor extraction/readback for reused meshIDs.
    // -----------------------------------------------------------------
    for (uint32_t instanceID : instanceIDs)
    {
        const auto& instance =
            mpScene->getGeometryInstance(instanceID);

        const MeshID meshID(instance.geometryID);
        const uint32_t meshIDValue = meshID.get();

        if (instance.globalMatrixID >= globalMatrices.size())
        {
            xatlas::Destroy(atlas);

            FALCOR_THROW(
                "Invalid globalMatrixID={} for instanceID={}.",
                instance.globalMatrixID,
                instanceID
            );
        }

        const auto cacheIt =
            mMeshGeometryCache.find(meshIDValue);

        if (cacheIt == mMeshGeometryCache.end())
        {
            xatlas::Destroy(atlas);
            FALCOR_THROW(
                "Missing cached mesh geometry for meshID={}.",
                meshIDValue
            );
        }

        const CachedMeshGeometry& geometry =
            cacheIt->second;

        // -------------------------------------------------------------
        // Transform this particular INSTANCE to world space before giving
        // its geometry to xatlas.
        //
        // geometry.positions is mesh-local.
        // globalMatrices[instance.globalMatrixID] is the same world
        // transform Falcor associates with this GeometryInstanceData.
        // -------------------------------------------------------------
        const float4x4& worldMat =
            globalMatrices[instance.globalMatrixID];

        std::vector<float3> worldPositions(
            geometry.positions.size()
        );

        for (size_t vertexIndex = 0;
            vertexIndex < geometry.positions.size();
            ++vertexIndex)
        {
            const float4 p =
                mul(
                    worldMat,
                    float4(
                        geometry.positions[vertexIndex],
                        1.0f
                    )
                );

            worldPositions[vertexIndex] =
                p.xyz();
        }

        // -------------------------------------------------------------
        // Same topology/indices as before, only geometry metric supplied
        // to xatlas is now the actual transformed instance geometry.
        // -------------------------------------------------------------
        xatlas::MeshDecl meshDecl{};

        meshDecl.vertexCount =
            uint32_t(worldPositions.size());

        meshDecl.vertexPositionData =
            worldPositions.data();

        meshDecl.vertexPositionStride =
            sizeof(float3);

        meshDecl.indexCount =
            uint32_t(geometry.indices.size());

        meshDecl.indexData =
            geometry.indices.data();

        meshDecl.indexFormat =
            xatlas::IndexFormat::UInt32;

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
                meshIDValue,
                errorString
            );
        }

        InputMeshInfo input;
        input.instanceID = instanceID;
        input.meshID = meshIDValue;
        input.triangleCount =
            uint32_t(geometry.triangles.size());

        submittedTriangles += input.triangleCount;
        inputMeshes.push_back(input);
    }

    logInfo(
        "xatlas global input submitted: meshes={} triangles={}. "
        "Waiting for AddMesh work...",
        inputMeshes.size(),
        submittedTriangles
    );

    xatlas::AddMeshJoin(atlas);

    logInfo(
        "xatlas global AddMeshJoin complete. Starting ComputeCharts()..."
    );

    xatlas::ChartOptions chartOptions{};
    xatlas::ComputeCharts(atlas, chartOptions);

    logInfo(
        "xatlas global ComputeCharts complete. Starting PackCharts() "
        "with targetResolution={} and automatic texel density...",
        targetResolution
    );

    // -----------------------------------------------------------------
    // Important xatlas mode.
    //
    // xatlas documentation:
    // - texelsPerUnit == 0 asks xatlas to estimate the density.
    // - resolution != 0 provides the resolution that estimate should
    //   approximately match.
    //
    // This is fundamentally different from our old fixed-density mode
    // (resolution > 0 AND texelsPerUnit > 0), which explicitly allowed
    // xatlas to create many fixed-resolution sub-atlases.
    // -----------------------------------------------------------------
    xatlas::PackOptions packOptions{};
    packOptions.resolution = targetResolution;
    packOptions.texelsPerUnit = 0.0f;

    // Keep lightmap-safe sampling margins.
    packOptions.padding = 4;
    packOptions.bilinear = true;

    // Officially documented to reduce the number of possible chart
    // locations and improve packing speed.
    packOptions.blockAlign = true;

    // Random packing is much faster than brute-force packing.
    packOptions.bruteForce = false;

    packOptions.createImage = false;
    packOptions.rotateChartsToAxis = true;
    packOptions.rotateCharts = true;

    xatlas::PackCharts(atlas, packOptions);

    if (atlas->meshCount != inputMeshes.size())
    {
        const uint32_t outputMeshCount = atlas->meshCount;
        const size_t expectedMeshCount = inputMeshes.size();

        xatlas::Destroy(atlas);

        FALCOR_THROW(
            "xatlas mesh count mismatch after PackCharts(). "
            "Expected {}, got {}.",
            expectedMeshCount,
            outputMeshCount
        );
    }

    if (atlas->atlasCount == 0 ||
        atlas->width == 0 ||
        atlas->height == 0)
    {
        xatlas::Destroy(atlas);
        FALCOR_THROW("xatlas generated no valid atlas pages.");
    }

    outChosenTexelsPerUnit = atlas->texelsPerUnit;

    logInfo(
        "xatlas global result: meshes={} charts={} pages={} "
        "size={}x{} chosenTexelsPerUnit={}.",
        atlas->meshCount,
        atlas->chartCount,
        atlas->atlasCount,
        atlas->width,
        atlas->height,
        atlas->texelsPerUnit
    );

    for (uint32_t pageIndex = 0;
        pageIndex < atlas->atlasCount;
        ++pageIndex)
    {
        logInfo(
            "xatlas page {} utilization={:.2f}%.",
            pageIndex,
            atlas->utilization[pageIndex] * 100.f
        );
    }

    std::vector<AtlasPageData> pages(atlas->atlasCount);

    for (uint32_t pageIndex = 0;
        pageIndex < atlas->atlasCount;
        ++pageIndex)
    {
        AtlasPageData& page = pages[pageIndex];

        page.pageIndex = pageIndex;
        page.width = atlas->width;
        page.height = atlas->height;
        page.outputPath =
            "Bistro_AtlasPage_" +
            std::to_string(pageIndex) +
            ".exr";
        //page.outputPath =
        //    "Room_AtlasPage_" +
        //    std::to_string(pageIndex) +
        //    ".exr";
    }

    uint64_t packedTriangleCount = 0;
    uint64_t unatlasedTriangleCount = 0;

    // -----------------------------------------------------------------
    // Preserve the face mapping already validated in earlier tests:
    // source face f -> output indices f*3 + corner -> xatlas Vertex.
    // Vertex::xref must map back to the original source vertex.
    // -----------------------------------------------------------------
    for (uint32_t meshIndex = 0;
        meshIndex < atlas->meshCount;
        ++meshIndex)
    {
        const xatlas::Mesh& xaMesh =
            atlas->meshes[meshIndex];

        const InputMeshInfo& input =
            inputMeshes[meshIndex];

        const auto cacheIt =
            mMeshGeometryCache.find(input.meshID);

        if (cacheIt == mMeshGeometryCache.end())
        {
            xatlas::Destroy(atlas);
            FALCOR_THROW(
                "Cached source geometry disappeared for meshID={}.",
                input.meshID
            );
        }

        const CachedMeshGeometry& geometry =
            cacheIt->second;

        const uint32_t triangleCount =
            input.triangleCount;

        if (xaMesh.indexCount != triangleCount * 3)
        {
            const uint32_t outputIndexCount =
                xaMesh.indexCount;

            xatlas::Destroy(atlas);

            FALCOR_THROW(
                "xatlas index-count mismatch for instanceID={}. "
                "Expected {}, got {}.",
                input.instanceID,
                triangleCount * 3,
                outputIndexCount
            );
        }

        for (uint32_t triangleID = 0;
            triangleID < triangleCount;
            ++triangleID)
        {
            const uint32_t firstIndex =
                triangleID * 3;

            const xatlas::Vertex& v0 =
                xaMesh.vertexArray[
                    xaMesh.indexArray[firstIndex + 0]
                ];
            const xatlas::Vertex& v1 =
                xaMesh.vertexArray[
                    xaMesh.indexArray[firstIndex + 1]
                ];
            const xatlas::Vertex& v2 =
                xaMesh.vertexArray[
                    xaMesh.indexArray[firstIndex + 2]
                ];

            const uint3& sourceTriangle =
                geometry.triangles[triangleID];

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

            const int32_t atlasIndex =
                v0.atlasIndex;

            if (v1.atlasIndex != atlasIndex ||
                v2.atlasIndex != atlasIndex)
            {
                xatlas::Destroy(atlas);

                FALCOR_THROW(
                    "xatlas produced a cross-page face for "
                    "instanceID={} triangleID={}: ({}, {}, {}).",
                    input.instanceID,
                    triangleID,
                    v0.atlasIndex,
                    v1.atlasIndex,
                    v2.atlasIndex
                );
            }

            // Same behavior as the official examples: skip faces that xatlas
            // could not atlas, such as degenerate geometry.
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
                pages[uint32_t(atlasIndex)];

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

            AtlasPageTriangleData triangleData;
            triangleData.instanceID =
                input.instanceID;
            triangleData.meshID =
                input.meshID;
            triangleData.triangleID =
                triangleID;
            triangleData.uv =
                uv;

            page.triangles.push_back(
                std::move(triangleData)
            );

            ++packedTriangleCount;
        }
    }

    xatlas::Destroy(atlas);

    uint64_t storedTriangleCount = 0;

    for (const AtlasPageData& page : pages)
    {
        storedTriangleCount +=
            page.triangles.size();

        logInfo(
            "Atlas page {} stored: size={}x{} triangles={} output='{}'.",
            page.pageIndex,
            page.width,
            page.height,
            page.triangles.size(),
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
        "Global atlas conversion complete: packedTriangles={} "
        "unatlasedTriangles={} pages={}.",
        packedTriangleCount,
        unatlasedTriangleCount,
        pages.size()
    );

    return pages;
}

void BakeLightMapXAtlas::buildAndSaveAtlasMapping(
    const std::vector<uint32_t>& instanceIDs
)
{
    if (instanceIDs.empty())
        FALCOR_THROW("Cannot build atlas mapping from an empty instance list.");

    // Extract each unique Falcor mesh once. Reused scene instances share this
    // cached CPU geometry when being submitted to xatlas.
    buildMeshGeometryCache(instanceIDs);

    float chosenTexelsPerUnit = 0.f;

    // One global xatlas job, matching the official example structure:
    // all inputs -> ComputeCharts once -> PackCharts once.
    std::vector<AtlasPageData> pages =
        buildGlobalAtlas(
            instanceIDs,
            mAtlasResolution,
            chosenTexelsPerUnit
        );

    if (pages.empty())
        FALCOR_THROW("Global xatlas build produced no pages.");

    AtlasMappingHeader header;
    header.pageCount =
        uint32_t(pages.size());
    header.width =
        pages.front().width;
    header.height =
        pages.front().height;
    header.texelsPerUnit =
        chosenTexelsPerUnit;

    uint64_t totalRecordCount = 0;

    for (const AtlasPageData& page : pages)
    {
        if (page.width != header.width ||
            page.height != header.height)
        {
            FALCOR_THROW(
                "Atlas page dimension mismatch. Expected {}x{}, "
                "got {}x{} on page {}.",
                header.width,
                header.height,
                page.width,
                page.height,
                page.pageIndex
            );
        }

        totalRecordCount +=
            page.triangles.size();
    }

    header.triangleRecordCount =
        totalRecordCount;

    // Write only after xatlas succeeds. This avoids leaving a partially valid
    // mapping header if charting/packing fails.
    std::ofstream mappingFile(
        kAtlasMappingFile,
        std::ios::binary | std::ios::trunc
    );

    if (!mappingFile)
        FALCOR_THROW(
            "Failed to open '{}' for writing.",
            kAtlasMappingFile
        );

    mappingFile.write(
        reinterpret_cast<const char*>(&header),
        sizeof(header)
    );

    for (const AtlasPageData& page : pages)
    {
        for (const AtlasPageTriangleData& triangle : page.triangles)
        {
            const TriangleLightmapUV& uv =
                triangle.uv;

            AtlasMappingRecord record;
            record.pageIndex =
                page.pageIndex;
            record.instanceID =
                triangle.instanceID;
            record.meshID =
                triangle.meshID;
            record.triangleID =
                triangle.triangleID;

            record.uv[0] = uv.uv0.x;
            record.uv[1] = uv.uv0.y;
            record.uv[2] = uv.uv1.x;
            record.uv[3] = uv.uv1.y;
            record.uv[4] = uv.uv2.x;
            record.uv[5] = uv.uv2.y;

            mappingFile.write(
                reinterpret_cast<const char*>(&record),
                sizeof(record)
            );
        }
    }

    if (!mappingFile)
        FALCOR_THROW(
            "Failed while writing '{}'.",
            kAtlasMappingFile
        );

    mappingFile.close();

    saveAtlasManifest(
        header.pageCount,
        header.width,
        header.height,
        header.triangleRecordCount,
        header.texelsPerUnit
    );

    // xatlas is finished and the persistent mapping now owns everything needed
    // by the bake. Release the CPU mesh cache before allocating bake resources.
    mMeshGeometryCache.clear();

    logInfo(
        "Global atlas mapping complete: pages={} records={} size={}x{} "
        "chosenTexelsPerUnit={} mapping='{}'.",
        header.pageCount,
        header.triangleRecordCount,
        header.width,
        header.height,
        header.texelsPerUnit,
        kAtlasMappingFile
    );
}

std::vector<BakeLightMapXAtlas::AtlasPageData>
BakeLightMapXAtlas::loadAtlasPagesFromMapping() const
{
    std::ifstream mappingFile(
        kAtlasMappingFile,
        std::ios::binary
    );

    if (!mappingFile)
        FALCOR_THROW("Failed to open '{}' for reading.", kAtlasMappingFile);

    AtlasMappingHeader header;
    mappingFile.read(
        reinterpret_cast<char*>(&header),
        sizeof(header)
    );

    if (!mappingFile)
        FALCOR_THROW("Failed to read atlas mapping header.");

    static const char kExpectedMagic[8] =
    { 'X', 'A', 'L', 'M', 'A', 'P', '0', '1' };

    if (std::memcmp(header.magic, kExpectedMagic, 8) != 0)
        FALCOR_THROW("Invalid atlas mapping magic in '{}'.", kAtlasMappingFile);
    if (header.version != 1)
        FALCOR_THROW("Unsupported atlas mapping version {}.", header.version);
    if (header.pageCount == 0 || header.width == 0 || header.height == 0)
        FALCOR_THROW("Atlas mapping header contains invalid page dimensions.");

    std::vector<AtlasPageData> pages(header.pageCount);

    for (uint32_t pageIndex = 0;
        pageIndex < header.pageCount;
        ++pageIndex)
    {
        AtlasPageData& page = pages[pageIndex];
        page.pageIndex = pageIndex;
        page.width = header.width;
        page.height = header.height;
        page.outputPath =
            "Bistro_AtlasPage_" +
            std::to_string(pageIndex) +
            ".exr";
        //page.outputPath =
        //    "Room_AtlasPage_" +
        //    std::to_string(pageIndex) +
        //    ".exr";
    }

    for (uint64_t recordIndex = 0;
        recordIndex < header.triangleRecordCount;
        ++recordIndex)
    {
        AtlasMappingRecord record;
        mappingFile.read(
            reinterpret_cast<char*>(&record),
            sizeof(record)
        );

        if (!mappingFile)
        {
            FALCOR_THROW(
                "Atlas mapping ended early at record {} / {}.",
                recordIndex,
                header.triangleRecordCount
            );
        }

        if (record.pageIndex >= header.pageCount)
        {
            FALCOR_THROW(
                "Atlas mapping record {} has invalid pageIndex={}.",
                recordIndex,
                record.pageIndex
            );
        }

        AtlasPageTriangleData triangle;
        triangle.instanceID = record.instanceID;
        triangle.meshID = record.meshID;
        triangle.triangleID = record.triangleID;
        triangle.uv.uv0 = float2(record.uv[0], record.uv[1]);
        triangle.uv.uv1 = float2(record.uv[2], record.uv[3]);
        triangle.uv.uv2 = float2(record.uv[4], record.uv[5]);

        pages[record.pageIndex].triangles.push_back(
            std::move(triangle)
        );
    }

    logInfo(
        "Loaded persistent atlas mapping: pages={} records={} size={}x{} "
        "texelsPerUnit={}.",
        header.pageCount,
        header.triangleRecordCount,
        header.width,
        header.height,
        header.texelsPerUnit
    );

    return pages;
}

void BakeLightMapXAtlas::saveAtlasManifest(
    uint32_t pageCount,
    uint32_t width,
    uint32_t height,
    uint64_t triangleRecordCount,
    float texelsPerUnit
) const
{
    std::ofstream manifestFile(
        kAtlasManifestFile,
        std::ios::out | std::ios::trunc
    );

    if (!manifestFile)
        FALCOR_THROW("Failed to open '{}' for writing.", kAtlasManifestFile);

    manifestFile
        << "format=XALMAP01\n"
        << "version=1\n"
        << "mappingFile=" << kAtlasMappingFile << "\n"
        << "pageCount=" << pageCount << "\n"
        << "width=" << width << "\n"
        << "height=" << height << "\n"
        << "texelsPerUnit=" << texelsPerUnit << "\n"
        << "triangleRecordCount=" << triangleRecordCount << "\n";

    for (uint32_t pageIndex = 0;
        pageIndex < pageCount;
        ++pageIndex)
    {
        manifestFile
            << "page" << pageIndex
            << "=Bistro_AtlasPage_" << pageIndex
            << ".exr\n";
        //manifestFile
        //    << "page" << pageIndex
        //    << "=Room_AtlasPage_" << pageIndex
        //    << ".exr\n";
    }

    if (!manifestFile)
        FALCOR_THROW("Failed while writing '{}'.", kAtlasManifestFile);

    logInfo(
        "Saved atlas manifest '{}' for {} pages.",
        kAtlasManifestFile,
        pageCount
    );
}

void BakeLightMapXAtlas::createAtlasPageGpuBuffers(
    AtlasPageData& page
)
{
    const uint32_t triangleCount =
        uint32_t(page.triangles.size());

    if (triangleCount == 0)
        return;

    std::vector<uint32_t> instanceIDs(triangleCount);
    std::vector<uint32_t> triangleIDs(triangleCount);
    std::vector<TriangleLightmapUV> triangleUVs(triangleCount);

    for (uint32_t i = 0; i < triangleCount; ++i)
    {
        const AtlasPageTriangleData& triangle = page.triangles[i];
        instanceIDs[i] = triangle.instanceID;
        triangleIDs[i] = triangle.triangleID;
        triangleUVs[i] = triangle.uv;
    }

    page.pInstanceIDBuffer =
        mpDevice->createStructuredBuffer(
            sizeof(uint32_t),
            triangleCount,
            ResourceBindFlags::ShaderResource
        );
    page.pInstanceIDBuffer->setBlob(
        instanceIDs.data(),
        0,
        instanceIDs.size() * sizeof(uint32_t)
    );

    page.pTriangleIDBuffer =
        mpDevice->createStructuredBuffer(
            sizeof(uint32_t),
            triangleCount,
            ResourceBindFlags::ShaderResource
        );
    page.pTriangleIDBuffer->setBlob(
        triangleIDs.data(),
        0,
        triangleIDs.size() * sizeof(uint32_t)
    );

    page.pTriangleUVBuffer =
        mpDevice->createStructuredBuffer(
            sizeof(TriangleLightmapUV),
            triangleCount,
            ResourceBindFlags::ShaderResource
        );
    page.pTriangleUVBuffer->setBlob(
        triangleUVs.data(),
        0,
        triangleUVs.size() * sizeof(TriangleLightmapUV)
    );

    logInfo(
        "Created GPU raster buffers for atlas page {} ({} triangles).",
        page.pageIndex,
        triangleCount
    );
}

void BakeLightMapXAtlas::releaseAtlasPageGpuBuffers(
    AtlasPageData& page
)
{
    page.pInstanceIDBuffer = nullptr;
    page.pTriangleIDBuffer = nullptr;
    page.pTriangleUVBuffer = nullptr;
}

void BakeLightMapXAtlas::bakeSavedAtlasPages(
    RenderContext* pRenderContext
)
{
    std::vector<AtlasPageData> pages =
        loadAtlasPagesFromMapping();

    logInfo(
        "Starting lightmap bake from persistent mapping: pages={} "
        "samplesPerPage={}.",
        pages.size(),
        mBakeSampleCount
    );

    for (size_t pageIndex = 0;
        pageIndex < pages.size();
        ++pageIndex)
    {
        AtlasPageData& page = pages[pageIndex];

        logInfo(
            "Starting bake for atlas page {}/{} (globalPage={} triangles={})...",
            pageIndex + 1,
            pages.size(),
            page.pageIndex,
            page.triangles.size()
        );

        createAtlasPageGpuBuffers(page);
        bakeAtlasPage(pRenderContext, page);
        releaseAtlasPageGpuBuffers(page);

        // This page will never be needed again in this pass. Release its CPU
        // triangle mapping too, so memory falls as the bake advances.
        page.triangles.clear();
        page.triangles.shrink_to_fit();
    }

    logInfo(
        "Full lightmap bake complete: pages={} samplesPerPage={}.",
        pages.size(),
        mBakeSampleCount
    );
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

    const uint32_t rasterTriangleCount =
        uint32_t(page.triangles.size());

    // -----------------------------------------------------------------
    // Draw the entire physical page in one procedural draw. Each page-local
    // triangle fetches its source Falcor instance and triangle ID from the
    // parallel indirection buffers.
    // -----------------------------------------------------------------
    if (rasterTriangleCount > 0)
    {
        if (!page.pInstanceIDBuffer ||
            !page.pTriangleIDBuffer ||
            !page.pTriangleUVBuffer)
        {
            FALCOR_THROW(
                "Atlas page {} is missing GPU raster data.",
                page.pageIndex
            );
        }

        var["gInstanceIDs"] =
            page.pInstanceIDBuffer;
        var["gTriangleIDs"] =
            page.pTriangleIDBuffer;
        var["gTriangleUVs"] =
            page.pTriangleUVBuffer;

        const uint32_t vertexCount =
            rasterTriangleCount * 3;

        pRenderContext->draw(
            mpUVGraphicsState.get(),
            mpUVVars.get(),
            vertexCount,
            0
        );

        pPosTex->captureToFile(
            0,
            0,
            "DEBUG_XAtlas_Pos_bistro.exr",
            Bitmap::FileFormat::ExrFile,
            Bitmap::ExportFlags::Uncompressed
        );

        pNormTex->captureToFile(
            0,
            0,
            "DEBUG_XAtlas_Normal_bistro.exr",
            Bitmap::FileFormat::ExrFile,
            Bitmap::ExportFlags::Uncompressed
        );
    }

    logInfo(
        "Rasterized atlas page {}: size={}x{} triangles={} drawCalls=1.",
        page.pageIndex,
        mLightmapWidth,
        mLightmapHeight,
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
