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

const char kUVRasterFile[] =
"RenderPasses/BakeLightMapXAtlas/XAtlasUVRaster.slang";

const char kBakingFile[] =
"RenderPasses/BakeLightMapXAtlas/LightmapBakingSingle.rt.slang";

const char kExtractFile[] =
"RenderPasses/BakeLightMapXAtlas/ExtractTexelsSingle.cs.slang";

const char kNormalizeFile[] =
"RenderPasses/BakeLightMapXAtlas/NormalizeLightmapSingle.cs.slang";

extern "C" FALCOR_API_EXPORT void registerPlugin(Falcor::PluginRegistry& registry)
{
    registry.registerClass<RenderPass, BakeLightMapXAtlas>();
}

BakeLightMapXAtlas::BakeLightMapXAtlas(ref<Device> pDevice, const Properties& props) : RenderPass(pDevice)
{

}

Properties BakeLightMapXAtlas::getProperties() const
{
    return {};
}

RenderPassReflection BakeLightMapXAtlas::reflect(const CompileData& compileData)
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

void BakeLightMapXAtlas::execute(
    RenderContext* pRenderContext,
    const RenderData& renderData
)
{
    if (!mpScene)
    {
        return;
    }

    if (mCurrentTargetIndex >= mBakeTargets.size())
    {
        return;
    }

    const BakeTarget& target =
        mBakeTargets[mCurrentTargetIndex];

    // -------------------------------------------------------------
    // Prepare UV atlas, extract valid texels, allocate accumulation.
    // Runs only once per target.
    // -------------------------------------------------------------

    if (mNeedsPreparation)
    {
        prepareTarget(
            pRenderContext,
            target
        );

        if (mNumExtractedTexels == 0)
        {
            logWarning(
                "Target '{}' has no extracted texels. Skipping.",
                target.name
            );

            advanceTarget();
            return;
        }
    }

    // -------------------------------------------------------------
    // Accumulate one RT sample per frame.
    // -------------------------------------------------------------

    if (mCurrentSample < mTotalSamples)
    {
        traceOneSample(
            pRenderContext
        );

        return;
    }

    // -------------------------------------------------------------
    // All samples complete.
    // Normalize and save EXR once.
    // -------------------------------------------------------------

    if (mCurrentSample == mTotalSamples)
    {
        finalizeTarget(
            pRenderContext,
            target
        );

        logInfo(
            "Finished baking '{}' with {} samples.",
            target.name,
            mTotalSamples
        );

        advanceTarget();

        return;
    }
}

void BakeLightMapXAtlas::renderUI(Gui::Widgets& widget) {}
void BakeLightMapXAtlas::setScene(
    RenderContext* pRenderContext,
    const ref<Scene>& pScene
)
{
    resetBakingState();

    mpScene = pScene;

    if (!mpScene)
    {
        return;
    }

    createUVRasterProgram();
    createComputePasses();
    createRayTracingProgram(pRenderContext);

    buildBakeTargets();
    buildMeshAtlases(pRenderContext);
}

void BakeLightMapXAtlas::resetBakingState()
{
    mMeshAtlases.clear();
    mBakeTargets.clear();

    mCurrentTargetIndex = 0;
    mCurrentSample = 0;
    mNumExtractedTexels = 0;
    mNeedsPreparation = true;

    mpUVFbo = nullptr;
    mpTexelBuffer = nullptr;
    mpCounterBuffer = nullptr;
    mpAccumBuffer = nullptr;
    mpResultTex = nullptr;

    mpUVProgram = nullptr;
    mpUVVars = nullptr;
    mpUVGraphicsState = nullptr;

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
    rasterDesc.setFillMode(
        RasterizerState::FillMode::Solid
    );
    rasterDesc.setCullMode(
        RasterizerState::CullMode::None
    );

    mpRasterState =
        RasterizerState::create(rasterDesc);

    DepthStencilState::Desc depthDesc;
    depthDesc.setDepthEnabled(false);

    ref<DepthStencilState> pDepthState =
        DepthStencilState::create(depthDesc);

    mpUVGraphicsState =
        GraphicsState::create(mpDevice);

    mpUVGraphicsState->setProgram(mpUVProgram);
    mpUVGraphicsState->setRasterizerState(mpRasterState);
    mpUVGraphicsState->setDepthStencilState(pDepthState);
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
    rtDesc.addShaderModules(
        mpScene->getShaderModules()
    );
    rtDesc.addShaderLibrary(kBakingFile);
    rtDesc.setMaxTraceRecursionDepth(3);
    rtDesc.setMaxPayloadSize(128);
    rtDesc.setMaxAttributeSize(8);
    rtDesc.addTypeConformances(
        mpScene->getTypeConformances()
    );

    ref<RtBindingTable> sbt =
        RtBindingTable::create(
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
        mpScene->getGeometryIDs(
            Scene::GeometryType::TriangleMesh
        );

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
            FALCOR_THROW(
                "Unknown emissive light sampler type."
            );
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

void BakeLightMapXAtlas::buildBakeTargets()
{
    mBakeTargets.clear();
    //std::unordered_map<uint32_t, uint32_t> firstInstanceForMesh;

    //for (uint32_t instanceID = 0;
    //    instanceID < mpScene->getGeometryInstanceCount();
    //    ++instanceID)
    //{
    //    const auto& instance =
    //        mpScene->getGeometryInstance(instanceID);

    //    const uint32_t meshID =
    //        instance.geometryID;

    //    auto it =
    //        firstInstanceForMesh.find(meshID);

    //    if (it == firstInstanceForMesh.end())
    //    {
    //        firstInstanceForMesh[meshID] =
    //            instanceID;
    //    }
    //    else
    //    {
    //        logInfo(
    //            "Found repeated mesh: meshID={} instances={} and {}",
    //            meshID,
    //            it->second,
    //            instanceID
    //        );

    //        break;
    //    }
    //}
    const uint32_t testInstanceIDs[] =
    {
        1786,
        1787
    };

    for (uint32_t instanceID : testInstanceIDs)
    {
        if (instanceID >= mpScene->getGeometryInstanceCount())
        {
            logWarning(
                "Skipping invalid test instance ID {}.",
                instanceID
            );

            continue;
        }

        const auto& instance =
            mpScene->getGeometryInstance(instanceID);

        const MeshID meshID(
            instance.geometryID
        );

        const auto& mesh =
            mpScene->getMesh(meshID);

        BakeTarget target;

        target.name =
            "BistroTest_Instance" +
            std::to_string(instanceID);

        target.instanceID =
            instanceID;

        target.meshID =
            meshID;

        target.outputPath =
            target.name + ".exr";

        mBakeTargets.push_back(
            target
        );

        logInfo(
            "Added bake target '{}' "
            "instanceID={} meshID={} "
            "vertices={} triangles={}",
            target.name,
            target.instanceID,
            target.meshID,
            mesh.vertexCount,
            mesh.indexCount / 3
        );
    }
}

void BakeLightMapXAtlas::buildMeshAtlases(RenderContext* pRenderContext)
{
    mMeshAtlases.clear();

    for (const auto& target : mBakeTargets)
    {
        const uint32_t meshKey = target.meshID.get();

        // -------------------------------------------------------------
        // A mesh only needs one xatlas parameterization.
        // Multiple instances of the same mesh can reuse it.
        // -------------------------------------------------------------
        if (mMeshAtlases.find(meshKey) != mMeshAtlases.end())
        {
            logInfo(
                "Reusing xatlas meshID={} for instanceID={}",
                meshKey,
                target.instanceID
            );

            continue;
        }

        const auto& mesh = mpScene->getMesh(target.meshID);

        const uint32_t vertexCount = mesh.vertexCount;
        const uint32_t triangleCount = mesh.indexCount / 3;

        // =============================================================
        // STEP 1: Extract mesh data from Falcor
        // =============================================================

        ref<Buffer> pPositions =
            mpDevice->createStructuredBuffer(
                sizeof(float3),
                vertexCount,
                ResourceBindFlags::ShaderResource |
                ResourceBindFlags::UnorderedAccess
            );

        // Falcor requires the texcoord output buffer even though
        // xatlas does not need the original material UVs here.
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
            target.meshID,
            buffers
        );

        // =============================================================
        // STEP 2: Read extracted data back to CPU
        // =============================================================

        std::vector<float3> positions(vertexCount);
        std::vector<uint3> triangles(triangleCount);

        pPositions->getBlob(
            positions.data(),
            0,
            vertexCount * sizeof(float3)
        );

        pTriangleIndices->getBlob(
            triangles.data(),
            0,
            triangleCount * sizeof(uint3)
        );

        logInfo(
            "Extracted meshID={} vertices={} triangles={}",
            meshKey,
            vertexCount,
            triangleCount
        );

        // =============================================================
        // STEP 3: Flatten Falcor uint3 triangle indices
        //         into the array expected by xatlas
        // =============================================================

        std::vector<uint32_t> indices;
        indices.reserve(triangleCount * 3);

        for (const uint3& tri : triangles)
        {
            indices.push_back(tri.x);
            indices.push_back(tri.y);
            indices.push_back(tri.z);
        }

        // =============================================================
        // STEP 4: Create xatlas input
        // =============================================================

        xatlas::Atlas* atlas = xatlas::Create();

        xatlas::MeshDecl meshDecl{};

        meshDecl.vertexCount = vertexCount;

        meshDecl.vertexPositionData =
            positions.data();

        meshDecl.vertexPositionStride =
            sizeof(float3);

        meshDecl.indexCount =
            uint32_t(indices.size());

        meshDecl.indexData =
            indices.data();

        meshDecl.indexFormat =
            xatlas::IndexFormat::UInt32;

        // =============================================================
        // STEP 5: Add the mesh to xatlas
        // =============================================================

        const xatlas::AddMeshError error =
            xatlas::AddMesh(
                atlas,
                meshDecl
            );

        if (error != xatlas::AddMeshError::Success)
        {
            logError(
                "xatlas::AddMesh() failed for meshID={}: {}",
                meshKey,
                xatlas::StringForEnum(error)
            );

            xatlas::Destroy(atlas);

            FALCOR_THROW(
                "xatlas::AddMesh() failed."
            );
        }

        // =============================================================
        // STEP 6: Generate charts and pack them into the atlas
        // =============================================================

        xatlas::ChartOptions chartOptions{};

        xatlas::PackOptions packOptions{};

        // Temporary resolution for this first test.
        packOptions.resolution = 512;

        // Padding between UV charts.
        packOptions.padding = 4;

        // Reserve padding appropriate for bilinear filtering.
        packOptions.bilinear = true;

        // We don't need xatlas's debug image.
        packOptions.createImage = false;

        xatlas::Generate(
            atlas,
            chartOptions,
            packOptions
        );

        // =============================================================
        // STEP 7: Inspect xatlas output
        // =============================================================

        logInfo(
            "xatlas meshID={}: "
            "atlasCount={} "
            "size={}x{} "
            "charts={} "
            "meshes={}",
            meshKey,
            atlas->atlasCount,
            atlas->width,
            atlas->height,
            atlas->chartCount,
            atlas->meshCount
        );

        if (atlas->meshCount != 1)
        {
            const uint32_t outputMeshCount =
                atlas->meshCount;

            xatlas::Destroy(atlas);

            FALCOR_THROW(
                "Expected 1 xatlas output mesh, got {}.",
                outputMeshCount
            );
        }

        // Our current representation only stores one UV coordinate
        // per triangle corner. It does not yet store atlasIndex.
        if (atlas->atlasCount != 1)
        {
            const uint32_t outputAtlasCount =
                atlas->atlasCount;

            xatlas::Destroy(atlas);

            FALCOR_THROW(
                "Expected 1 xatlas atlas for this test, got {}.",
                outputAtlasCount
            );
        }

        const xatlas::Mesh& xaMesh =
            atlas->meshes[0];

        logInfo(
            "xatlas output: "
            "inputVertices={} "
            "outputVertices={} "
            "inputIndices={} "
            "outputIndices={} "
            "charts={}",
            vertexCount,
            xaMesh.vertexCount,
            uint32_t(indices.size()),
            xaMesh.indexCount,
            xaMesh.chartCount
        );

        // =============================================================
        // STEP 8: Print a few xatlas vertices
        // =============================================================

        const uint32_t vertexPrintCount =
            std::min(xaMesh.vertexCount, 10u);

        for (uint32_t i = 0;
            i < vertexPrintCount;
            ++i)
        {
            const xatlas::Vertex& v =
                xaMesh.vertexArray[i];

            logInfo(
                "  xaVertex[{}]: "
                "xref={} "
                "uv=({}, {}) "
                "atlas={} "
                "chart={}",
                i,
                v.xref,
                v.uv[0],
                v.uv[1],
                v.atlasIndex,
                v.chartIndex
            );
        }

        // =============================================================
        // STEP 8.5:
        //
        // Convert xatlas output into:
        //
        //      one TriangleLightmapUV per Falcor triangle
        //
        // At the same time, verify that xatlas preserved the triangle
        // and corner correspondence that we need later when indexing
        // this array using SV_PrimitiveID.
        // =============================================================

        if (xaMesh.indexCount != triangleCount * 3)
        {
            const uint32_t xaIndexCount =
                xaMesh.indexCount;

            xatlas::Destroy(atlas);

            FALCOR_THROW(
                "xatlas index count mismatch. Expected {}, got {}.",
                triangleCount * 3,
                xaIndexCount
            );
        }

        std::vector<TriangleLightmapUV> triangleUVs;
        triangleUVs.resize(triangleCount);

        const float invWidth =
            1.0f / float(atlas->width);

        const float invHeight =
            1.0f / float(atlas->height);

        bool mappingValid = true;

        for (uint32_t triangleID = 0;
            triangleID < triangleCount;
            ++triangleID)
        {
            // Original Falcor triangle.
            const uint3& originalTri =
                triangles[triangleID];

            // xatlas outputs an index buffer referencing its own
            // seam-split vertex array.
            const uint32_t xaIndex0 =
                xaMesh.indexArray[triangleID * 3 + 0];

            const uint32_t xaIndex1 =
                xaMesh.indexArray[triangleID * 3 + 1];

            const uint32_t xaIndex2 =
                xaMesh.indexArray[triangleID * 3 + 2];

            const xatlas::Vertex& xaV0 =
                xaMesh.vertexArray[xaIndex0];

            const xatlas::Vertex& xaV1 =
                xaMesh.vertexArray[xaIndex1];

            const xatlas::Vertex& xaV2 =
                xaMesh.vertexArray[xaIndex2];

            // ---------------------------------------------------------
            // xref refers back to the original Falcor vertex.
            //
            // Therefore this test verifies:
            //
            // Falcor triangle N:
            //      (v0, v1, v2)
            //
            // corresponds to xatlas triangle N:
            //      (xref0, xref1, xref2)
            //
            // with the same corner ordering.
            // ---------------------------------------------------------

            if (xaV0.xref != originalTri.x ||
                xaV1.xref != originalTri.y ||
                xaV2.xref != originalTri.z)
            {
                logError(
                    "xatlas triangle mapping mismatch at triangle {}: "
                    "Falcor=({}, {}, {}) "
                    "xatlas xref=({}, {}, {})",
                    triangleID,
                    originalTri.x,
                    originalTri.y,
                    originalTri.z,
                    xaV0.xref,
                    xaV1.xref,
                    xaV2.xref
                );

                mappingValid = false;
                break;
            }

            // ---------------------------------------------------------
            // xatlas UVs are in atlas texel coordinates.
            //
            // Convert them to normalized [0,1] coordinates.
            // ---------------------------------------------------------

            TriangleLightmapUV triUV;

            triUV.uv0 = float2(
                xaV0.uv[0] * invWidth,
                xaV0.uv[1] * invHeight
            );

            triUV.uv1 = float2(
                xaV1.uv[0] * invWidth,
                xaV1.uv[1] * invHeight
            );

            triUV.uv2 = float2(
                xaV2.uv[0] * invWidth,
                xaV2.uv[1] * invHeight
            );

            triangleUVs[triangleID] =
                triUV;
        }

        if (!mappingValid)
        {
            xatlas::Destroy(atlas);

            FALCOR_THROW(
                "xatlas triangle ordering does not match "
                "Falcor mesh ordering."
            );
        }

        logInfo(
            "xatlas triangle mapping verified for {} triangles.",
            triangleCount
        );

        // =============================================================
        // STEP 8.6: Print first few triangle UV triplets
        // =============================================================

        const uint32_t triangleUVPrintCount =
            std::min(triangleCount, 5u);

        for (uint32_t i = 0;
            i < triangleUVPrintCount;
            ++i)
        {
            const TriangleLightmapUV& uv =
                triangleUVs[i];

            logInfo(
                "  triangleUV[{}]: "
                "({}, {}) "
                "({}, {}) "
                "({}, {})",
                i,
                uv.uv0.x,
                uv.uv0.y,
                uv.uv1.x,
                uv.uv1.y,
                uv.uv2.x,
                uv.uv2.y
            );
        }

 // =============================================================
// STEP 9: Create GPU triangle-UV buffer and store atlas data
// =============================================================

        MeshAtlasData data;

        data.meshID =
            target.meshID;

        data.triangleCount =
            triangleCount;

        data.width =
            atlas->width;

        data.height =
            atlas->height;

        // -------------------------------------------------------------
        // Create one StructuredBuffer element per triangle.
        //
        // TriangleLightmapUV contains:
        //
        //     float2 uv0;
        //     float2 uv1;
        //     float2 uv2;
        //
        // so triangle N can later be accessed directly using
        // SV_PrimitiveID.
        // -------------------------------------------------------------

        data.pTriangleUVBuffer =
            mpDevice->createStructuredBuffer(
                sizeof(TriangleLightmapUV),
                triangleCount,
                ResourceBindFlags::ShaderResource
            );

        data.pTriangleUVBuffer->setBlob(
            triangleUVs.data(),
            0,
            triangleUVs.size() * sizeof(TriangleLightmapUV)
        );

        // Keep CPU copy as well for now.
        // This is useful while we're debugging the implementation.
        data.triangleUVs =
            std::move(triangleUVs);

        logInfo(
            "Created triangle UV buffer for meshID={} "
            "triangles={} elementSize={} totalBytes={}",
            meshKey,
            triangleCount,
            sizeof(TriangleLightmapUV),
            triangleCount * sizeof(TriangleLightmapUV)
        );

        mMeshAtlases.emplace(
            meshKey,
            std::move(data)
        );

        // =============================================================
        // STEP 10: Done with xatlas CPU objects
        // =============================================================

        xatlas::Destroy(atlas);
    }
    // =============================================================
    // STEP 11: Assign atlas dimensions to each bake target
    // =============================================================
    for (auto& target : mBakeTargets)
    {
        const auto atlasIt =
            mMeshAtlases.find(target.meshID.get());

        if (atlasIt == mMeshAtlases.end())
        {
            FALCOR_THROW(
                "No mesh atlas found for bake target '{}' (meshID={}).",
                target.name,
                target.meshID.get()
            );
        }

        const MeshAtlasData& atlasData =
            atlasIt->second;

        target.width =
            atlasData.width;

        target.height =
            atlasData.height;

        logInfo(
            "Bake target '{}' uses lightmap size {}x{}",
            target.name,
            target.width,
            target.height
        );
    }
}

void BakeLightMapXAtlas::prepareTarget(
    RenderContext* pRenderContext,
    const BakeTarget& target
)
{
    // -------------------------------------------------------------
    // Find the xatlas data belonging to this target's mesh.
    // -------------------------------------------------------------

    auto atlasIt =
        mMeshAtlases.find(target.meshID.get());

    if (atlasIt == mMeshAtlases.end())
    {
        FALCOR_THROW(
            "No MeshAtlasData found for target '{}'.",
            target.name
        );
    }

    const MeshAtlasData& atlasData =
        atlasIt->second;

    if (!atlasData.pTriangleUVBuffer)
    {
        FALCOR_THROW(
            "Target '{}' has no triangle UV buffer.",
            target.name
        );
    }

    // -------------------------------------------------------------
    // Lightmap resolution comes directly from xatlas.
    // For our current mesh this is 118 x 115.
    // -------------------------------------------------------------

    mLightmapWidth =
        target.width;

    mLightmapHeight =
        target.height;

    if (mLightmapWidth == 0 ||
        mLightmapHeight == 0)
    {
        FALCOR_THROW(
            "Invalid lightmap size {}x{} for target '{}'.",
            mLightmapWidth,
            mLightmapHeight,
            target.name
        );
    }

    // -------------------------------------------------------------
    // Create textures that store world-space position and normal
    // for every valid lightmap texel.
    // -------------------------------------------------------------

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

    pPosTex->setName("XAtlas_PosTex");

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

    pNormTex->setName("XAtlas_NormTex");

    // -------------------------------------------------------------
    // Create UV-space FBO.
    //
    // PS output:
    //   SV_Target0 -> world position
    //   SV_Target1 -> world normal
    // -------------------------------------------------------------

    mpUVFbo =
        Fbo::create(mpDevice);

    mpUVFbo->attachColorTarget(
        pPosTex,
        0
    );

    mpUVFbo->attachColorTarget(
        pNormTex,
        1
    );

    // -------------------------------------------------------------
    // Set rasterization viewport to actual xatlas dimensions.
    // -------------------------------------------------------------

    GraphicsState::Viewport viewport(
        0.f,
        0.f,
        float(mLightmapWidth),
        float(mLightmapHeight),
        0.f,
        1.f
    );

    mpUVGraphicsState->setViewport(
        0,
        viewport
    );

    mpUVGraphicsState->setFbo(
        mpUVFbo
    );

    // Start from empty textures.
    pRenderContext->clearFbo(
        mpUVFbo.get(),
        float4(0.f),
        1.f,
        0,
        FboAttachmentType::Color
    );

    // -------------------------------------------------------------
    // Bind xatlas data.
    // -------------------------------------------------------------

    auto var =
        mpUVVars->getRootVar();

    mpScene->bindShaderData(
        var["gScene"]
    );

    var["gTriangleUVs"] =
        atlasData.pTriangleUVBuffer;

    var["BakeCB"]["gReceiverInstanceID"] =
        target.instanceID;

    var["BakeCB"]["gTriangleCount"] =
        atlasData.triangleCount;

    ref<Vao> pVao =
        mpScene->getMeshVao();

    if (!pVao)
    {
        pVao =
            mpScene->getMeshVao16();
    }

    if (!pVao)
    {
        FALCOR_THROW(
            "Scene has no mesh VAO available for UV rasterization."
        );
    }

    mpUVGraphicsState->setVao(
        pVao
    );

    const uint32_t drawVertexCount =
        atlasData.triangleCount * 3;

    pRenderContext->draw(
        mpUVGraphicsState.get(),
        mpUVVars.get(),
        drawVertexCount,
        0
    );

    //var["gTriangleUVs"] =
    //    atlasData.pTriangleUVBuffer;

    //var["BakeCB"]["gReceiverInstanceID"] =
    //    target.instanceID;

    //var["BakeCB"]["gTriangleCount"] =
    //    atlasData.triangleCount;

    // -------------------------------------------------------------
    // Rasterize.
    //
    // Scene::rasterize() submits the complete scene, but gsMain()
    // discards every geometry instance except target.instanceID.
    // -------------------------------------------------------------

    //mpScene->rasterize(
    //    pRenderContext,
    //    mpUVGraphicsState.get(),
    //    mpUVVars.get(),
    //    mpRasterState,
    //    mpRasterState
    //);

    logInfo(
        "UV-rasterized target '{}' "
        "instanceID={} meshID={} size={}x{} "
        "triangles={} drawVertices={}",
        target.name,
        target.instanceID,
        target.meshID,
        mLightmapWidth,
        mLightmapHeight,
        atlasData.triangleCount,
        drawVertexCount
    );
    // =============================================================
// Validate UV rasterization by extracting occupied texels
// =============================================================

    const uint32_t totalTexels =
        mLightmapWidth * mLightmapHeight;

    // One possible TexelSample for every lightmap pixel.
    mpTexelBuffer =
        mpDevice->createStructuredBuffer(
            sizeof(TexelSample),
            totalTexels
        );

    mpTexelBuffer->setName("XAtlas_TexelBuffer");

    // Atomic counter used by ExtractTexelsSingle.cs.slang.
    mpCounterBuffer =
        mpDevice->createBuffer(
            sizeof(uint32_t),
            ResourceBindFlags::UnorderedAccess,
            MemoryType::DeviceLocal
        );

    mpCounterBuffer->setName("XAtlas_CounterBuffer");

    // Start counter at zero.
    pRenderContext->clearUAV(
        mpCounterBuffer->getUAV().get(),
        uint4(0)
    );

    // -------------------------------------------------------------
    // Bind extraction shader resources.
    // -------------------------------------------------------------

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
        uint2(
            mLightmapWidth,
            mLightmapHeight
        );

    // -------------------------------------------------------------
    // Run extraction.
    // -------------------------------------------------------------

    mpExtractPass->execute(
        pRenderContext,
        mLightmapWidth,
        mLightmapHeight
    );

    // -------------------------------------------------------------
    // Read number of covered texels.
    // -------------------------------------------------------------

    uint32_t extractedCount = 0;

    mpCounterBuffer->getBlob(
        &extractedCount,
        0,
        sizeof(uint32_t)
    );

    mNumExtractedTexels =
        extractedCount;

    logInfo(
        "Extracted {} / {} texels ({:.2f}% coverage) "
        "for target '{}'.",
        mNumExtractedTexels,
        totalTexels,
        100.f * float(mNumExtractedTexels) / float(totalTexels),
        target.name
    );

    // =============================================================
// Create irradiance accumulation buffer.
//
// IMPORTANT:
// TexelSample::texelIndex refers to the original atlas pixel,
// so this buffer must have one element per atlas texel,
// NOT one element per extracted texel.
// =============================================================


    mpAccumBuffer =
        mpDevice->createStructuredBuffer(
            sizeof(float4),
            totalTexels
        );

    mpAccumBuffer->setName(
        "XAtlas_IrradianceAccum"
    );

    // Start completely empty.
    pRenderContext->clearUAV(
        mpAccumBuffer->getUAV().get(),
        float4(0.f)
    );

    mCurrentSample = 0;

    logInfo(
        "Created accumulation buffer with {} texels.",
        totalTexels
    );

    mNeedsPreparation = false;
}

void BakeLightMapXAtlas::traceOneSample(
    RenderContext* pRenderContext
)
{
    if (!mpTexelBuffer)
    {
        FALCOR_THROW(
            "Texel buffer is null."
        );
    }

    if (!mpAccumBuffer)
    {
        FALCOR_THROW(
            "Accumulation buffer is null."
        );
    }

    if (mNumExtractedTexels == 0)
    {
        FALCOR_THROW(
            "No extracted texels available for ray tracing."
        );
    }

    auto rtVar =
        mpRtVars->getRootVar();

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

    const uint32_t threadCount =
        mNumExtractedTexels;

    mpScene->raytrace(
        pRenderContext,
        mpRtProgram.get(),
        mpRtVars,
        uint3(threadCount, 1, 1)
    );

    ++mCurrentSample;

    logInfo(
        "Completed RT sample {} for {} texels.",
        mCurrentSample,
        mNumExtractedTexels
    );
}

void BakeLightMapXAtlas::finalizeTarget(
    RenderContext* pRenderContext,
    const BakeTarget& target
)
{
    if (!mpAccumBuffer)
    {
        FALCOR_THROW(
            "Accumulation buffer is null."
        );
    }

    if (mLightmapWidth == 0 || mLightmapHeight == 0)
    {
        FALCOR_THROW(
            "Invalid lightmap size {}x{}.",
            mLightmapWidth,
            mLightmapHeight
        );
    }

    // -------------------------------------------------------------
    // Create normalized lightmap output texture.
    // -------------------------------------------------------------

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

    mpResultTex->setName("XAtlas_ResultTex");

    // Clear to black first.
    pRenderContext->clearUAV(
        mpResultTex->getUAV().get(),
        float4(0.f)
    );

    // -------------------------------------------------------------
    // Bind normalization resources.
    // -------------------------------------------------------------

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

    // -------------------------------------------------------------
    // Run normalization.
    // -------------------------------------------------------------

    mpNormalizePass->execute(
        pRenderContext,
        mLightmapWidth,
        mLightmapHeight
    );

    logInfo(
        "Normalized target '{}' into texture {}x{} using {} sample(s).",
        target.name,
        mLightmapWidth,
        mLightmapHeight,
        mCurrentSample
    );

    if (mpResultTex)
    {
        mpResultTex->captureToFile(
            0,
            0,
            target.outputPath,
            Bitmap::FileFormat::ExrFile,
            Bitmap::ExportFlags::Uncompressed
        );

        logInfo(
            "Saved test lightmap to '{}'.",
            target.outputPath
        );
    }
}

void BakeLightMapXAtlas::advanceTarget()
{
    ++mCurrentTargetIndex;

    mCurrentSample = 0;
    mNumExtractedTexels = 0;
    mNeedsPreparation = true;
}
