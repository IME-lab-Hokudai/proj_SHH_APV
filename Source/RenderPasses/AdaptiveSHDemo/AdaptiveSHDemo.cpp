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
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>
#include <vector>
#include "AdaptiveSHDemo.h"
#include "Rendering/Lights/EmissivePowerSampler.h"
#include "Rendering/Lights/EmissiveUniformSampler.h"
#include "Utils/Logger.h"
#include "../PrecomputeSHCoefficients/FDPlaneSettings.h"

#define PROBE_MODE_ADAPTIVE 0
#define PROBE_MODE_UNIFORM  1
 // CHANGE THIS LINE TO SWITCH MODES:
#define CURRENT_PROBE_MODE PROBE_MODE_UNIFORM
//#define CURRENT_PROBE_MODE PROBE_MODE_ADAPTIVE

#if CURRENT_PROBE_MODE == PROBE_MODE_ADAPTIVE
//const std::string loadFromFileName = "Seeded8DirectAbsErr5SubwayCorridorNoOpen.txt";
//const std::string loadFromFileName = "DirectAbsErr20DataScene.txt";
//const std::string loadFromFileName = "DirectAbsErr10DataScene.txt";
//const std::string loadFromFileName = "DirectAbsErr8p5N6DataScene.txt";
//const std::string loadFromFileName = "DirectAbsErr5DataScene.txt";
//const std::string loadFromFileName = "DirectAbsErr10DataSceneAvg.txt";
//const std::string loadFromFileName = "DirectAbsErr5DataSceneAvg.txt";
//const std::string loadFromFileName = "DirectAbsErr2DataSceneAvg.txt";
//const std::string loadFromFileName = "DirectAbsErr10ResidualScaleMetric.txt";
//const std::string loadFromFileName = "DirectAbsErr5ResidualScaleMetric.txt";
//const std::string loadFromFileName = "DirectAbsErr2ResidualScaleMetric.txt";
//const std::string loadFromFileName = "DirectAbsErr1DataScene.txt";
//const std::string loadFromFileName = "DirectAbsErr7p95DataSceneIrrMetric.txt";
//const std::string loadFromFileName = "DirectAbsErr4p25DataSceneIrrMetric.txt";
//const std::string loadFromFileName = "DirectAbsErr2p1DataScenePresentativeNormalMetricAvg.txt";

//const std::string loadFromFileName = "DirectAbsErr2ResidualScaleV2MetricCornellThinSlab.txt"; //note this is actually shadow boundary scene
//const std::string loadFromFileName = "DirectAbsErr2N5HessianMetricThinSlabScene4096spp.txt";
//const std::string loadFromFileName = "DirectAbsErr2N5EdgeGradientMetricThinSlabScene4096spp.txt";
//const std::string loadFromFileName = "DirectAbsErr2HessianMetricCornellThinSlab.txt";
//const std::string loadFromFileName = "DirectAbsErr2N6EdgeGradientMetricDataScene4096spp.txt";
//const std::string loadFromFileName = "DirectAbsErr2EdgeMetricCornellThinSlabV2.txt";
//const std::string loadFromFileName = "DirectAbsErr2HessianMetricCornellThinSlabV2.txt";
//const std::string loadFromFileName = "DirectAbsErr2N6HessianMetricDataScene4096spp.txt";

//const std::string loadFromFileName = "BistroInterior_Wine_BF_HessianRel_S8_memU32.txt";
const std::string loadFromFileName = "SceneForData_BF_HessianRel_S8_memU32.txt";
//const std::string loadFromFileName = "DirectAbsErr100HessianMetricBistroNew.txt";

//const std::string loadFromFileName = "U64BistroCubeGrid.txt";
//const std::string loadFromFileName = "DirectAbsErr100EGCMetricBistroNew.txt";
//const std::string loadFromFileName = "DirectAbsErr50HessianMetricBistro.txt";
//const std::string loadFromFileName = "DirectAbsErr100HessianMetricCubeBistro.txt";
//const std::string loadFromFileName = "DirectAbsErr100EGCMetricCubeBistro.txt";
//const std::string loadFromFileName = "DirectAbsErr2EGCShadowBoundary.txt";
//const std::string loadFromFileName = "DirectAbsErr2HessianShadowBoundary.txt";
const char kShaderFile[] = "RenderPasses/AdaptiveSHDemo/AdaptiveGridShaderXAtlasTest.slang";
//const char kShaderFile[] = "RenderPasses/AdaptiveSHDemo/AdaptiveGridShader.slang";

#else
//const std::string loadFromFileName = "U32DataScene.txt";
//const std::string loadFromFileName = "U32DataScene_4096spp_4spt.txt";
//const std::string loadFromFileName = "U64CornellShadowBoundaryScene.txt";
//const std::string loadFromFileName = "U32CornellShadowBoundaryScene.txt";
//const std::string loadFromFileName = "U64DataScene.txt";
//const std::string loadFromFileName = "U64DataScene_4096spp.txt";
//const std::string loadFromFileName = "U64Bistro.txt";
//const std::string loadFromFileName = "U64BistroCubeGrid.txt"; 
//const std::string loadFromFileName = "U64Shadowboundary.txt"; 
//const std::string loadFromFileName = "BistroInterior_Wine_U32.txt"; 
//const std::string loadFromFileName = "BistroInterior_Wine_E10_U32.txt"; // emissive x10 test
const std::string loadFromFileName = "SceneForData_U32.txt"; 
//const std::string loadFromFileName = "BistroInterior_Wine_E10_U32.txt"; 
//const std::string loadFromFileName = "BistroInterior_Wine_E100_U32.txt"; 
//const std::string loadFromFileName = "U8DataScene.txt";
//const std::string loadFromFileName = "U64BistroCubeGridKPathAllDiffuse.txt";
const char kShaderFile[] = "RenderPasses/AdaptiveSHDemo/UniformGridShader.slang";
#endif

const char kDynamicPassFile[] = "RenderPasses/AdaptiveSHDemo/DynamicPass.slang";
const char kDynamicFilterFile[] = "RenderPasses/AdaptiveSHDemo/DynamicFilter.cs.slang";
const char kCompositePassFile[] = "RenderPasses/AdaptiveSHDemo/CompositeDynamic.slang";
extern "C" FALCOR_API_EXPORT void registerPlugin(Falcor::PluginRegistry& registry)
{
    registry.registerClass<RenderPass, AdaptiveSHDemo>();
}

AdaptiveSHDemo::AdaptiveSHDemo(ref<Device> pDevice, const Properties& props) : RenderPass(pDevice)
{
    mEnableGlass = props.get<bool>("runtimeGlass", mEnableGlass);
    mShInterpolation = std::min(props.get<uint32_t>("shInterpolation", mShInterpolation), 2u);
    mShSeam = std::min(props.get<uint32_t>("shSeam", mShSeam), 2u);
    mShNonNegative = props.get<bool>("shNonNegative", mShNonNegative);
    mGlassRefractionStrength = std::clamp(props.get<float>("glassRefraction", mGlassRefractionStrength), 0.f, 3.f);
    mLiquidOpticalDepth = std::clamp(props.get<float>("liquidOpticalDepth", mLiquidOpticalDepth), 0.f, 10.0f);
    mpFbo = Fbo::create(mpDevice);
    Sampler::Desc samplerDesc;
    samplerDesc.setFilterMode(TextureFilteringMode::Linear, TextureFilteringMode::Linear, TextureFilteringMode::Linear);
    mpLinearSampler = mpDevice->createSampler(samplerDesc);
}

namespace
{
    constexpr float kPi = 3.14159265358979323846f;

    //constexpr const char* kBistroAtlasMappingFile =
    //    "Bistro_AtlasMapping.bin";
    //constexpr const char* kBistroAtlasManifestFile =
    //    "Bistro_AtlasManifest.txt";
    //
    constexpr const char* kBistroAtlasMappingFile =
        "Room_AtlasMapping.bin";
    constexpr const char* kBistroAtlasManifestFile =
        "Room_AtlasManifest.txt";
    
    
    //constexpr const char* kBistroAtlasMappingFile =
    //    "ShadowBoundary_AtlasMapping.bin";
    //constexpr const char* kBistroAtlasManifestFile =
    //    "ShadowBoundary_AtlasManifest.txt";
    constexpr uint32_t kBistroTestMaxPages = 10;

#pragma pack(push, 1)
    struct BistroAtlasMappingHeader
    {
        char magic[8];
        uint32_t version;
        uint32_t pageCount;
        uint32_t width;
        uint32_t height;
        float texelsPerUnit;
        uint64_t triangleRecordCount;
    };

    struct BistroAtlasMappingRecord
    {
        uint32_t pageIndex;
        uint32_t instanceID;
        uint32_t meshID;
        uint32_t triangleID;
        float uv[6];
    };
#pragma pack(pop)

    static_assert(sizeof(BistroAtlasMappingHeader) == 36);
    static_assert(sizeof(BistroAtlasMappingRecord) == 40);

    std::array<ref<Texture>, kBistroTestMaxPages> spBistroAtlasPages;

    ref<Buffer> spBistroInstanceRanges;
    ref<Buffer> spBistroPageIndices;
    ref<Buffer> spBistroUV01;
    ref<Buffer> spBistroUV2;

    uint32_t sBistroPageCount = 0;
    uint32_t sBistroMaxInstanceID = 0;
    bool sBistroLightmapEnabled = false;

    float degreesToRadians(float degrees)
    {
        return degrees * kPi / 180.0f;
    }

    float radiansToDegrees(float radians)
    {
        return radians * 180.0f / kPi;
    }

    std::vector<std::string> readBistroAtlasPageFiles(
        uint32_t pageCount
    )
    {
        std::vector<std::string> pageFiles(pageCount);

        std::ifstream manifest(kBistroAtlasManifestFile);

        if (manifest)
        {
            std::string line;

            while (std::getline(manifest, line))
            {
                // We only want entries of the form:
                //
                //     page0=Bistro_AtlasPage_0.exr
                //     page1=Bistro_AtlasPage_1.exr
                //
                // Do NOT interpret "pageCount=..." as a page entry.
                if (line.size() < 6)
                    continue;

                if (line.rfind("page", 0) != 0)
                    continue;

                const char firstIndexChar = line[4];

                if (firstIndexChar < '0' ||
                    firstIndexChar > '9')
                {
                    continue;
                }

                const size_t equalsPos =
                    line.find('=');

                if (equalsPos == std::string::npos ||
                    equalsPos <= 4)
                {
                    continue;
                }

                const std::string indexString =
                    line.substr(
                        4,
                        equalsPos - 4
                    );

                uint32_t pageIndex = 0;

                try
                {
                    pageIndex =
                        uint32_t(std::stoul(indexString));
                }
                catch (...)
                {
                    logWarning(
                        "Ignoring malformed atlas manifest entry: '{}'.",
                        line
                    );

                    continue;
                }

                if (pageIndex >= pageCount)
                    continue;

                pageFiles[pageIndex] =
                    line.substr(equalsPos + 1);
            }
        }

        // Fall back to the baker's normal file naming scheme.
        for (uint32_t pageIndex = 0;
            pageIndex < pageCount;
            ++pageIndex)
        {
            if (pageFiles[pageIndex].empty())
            {
                pageFiles[pageIndex] =
                    "Bistro_AtlasPage_" +
                    std::to_string(pageIndex) +
                    ".exr";
                //pageFiles[pageIndex] =
                //    "Room_AtlasPage_" +
                //    std::to_string(pageIndex) +
                //    ".exr";
            }
        }

        return pageFiles;
    }

    void loadBistroXAtlasLightmapTest(
        const ref<Device>& pDevice,
        const ref<Scene>& pScene
    )
    {
        sBistroLightmapEnabled = false;

        for (auto& page : spBistroAtlasPages) page = nullptr;
        spBistroInstanceRanges = nullptr;
        spBistroPageIndices = nullptr;
        spBistroUV01 = nullptr;
        spBistroUV2 = nullptr;

        std::ifstream mappingFile(
            kBistroAtlasMappingFile,
            std::ios::binary
        );

        if (!mappingFile)
        {
            logWarning(
                "Bistro xatlas test disabled: cannot open '{}'.",
                kBistroAtlasMappingFile
            );
            return;
        }

        BistroAtlasMappingHeader header{};

        mappingFile.read(
            reinterpret_cast<char*>(&header),
            sizeof(header)
        );

        const char expectedMagic[8] =
        {
            'X', 'A', 'L', 'M', 'A', 'P', '0', '1'
        };

        if (!mappingFile ||
            std::memcmp(
                header.magic,
                expectedMagic,
                sizeof(expectedMagic)
            ) != 0)
        {
            FALCOR_THROW(
                "Invalid Bistro atlas mapping header in '{}'.",
                kBistroAtlasMappingFile
            );
        }

        if (header.version != 1)
        {
            FALCOR_THROW(
                "Unsupported Bistro atlas mapping version {}.",
                header.version
            );
        }

        if (header.pageCount == 0 ||
            header.pageCount > kBistroTestMaxPages)
        {
            FALCOR_THROW(
                "Bistro xatlas TEST shader supports 1-{} pages, "
                "but mapping contains {} page(s).",
                kBistroTestMaxPages,
                header.pageCount
            );
        }

        if (header.triangleRecordCount >
            uint64_t(std::numeric_limits<size_t>::max()))
        {
            FALCOR_THROW(
                "Bistro atlas mapping record count is too large."
            );
        }

        std::vector<BistroAtlasMappingRecord> records(
            size_t(header.triangleRecordCount)
        );

        if (!records.empty())
        {
            mappingFile.read(
                reinterpret_cast<char*>(records.data()),
                std::streamsize(
                    records.size() *
                    sizeof(BistroAtlasMappingRecord)
                )
            );
        }

        if (!mappingFile)
        {
            FALCOR_THROW(
                "Failed while reading '{}'.",
                kBistroAtlasMappingFile
            );
        }

        uint32_t maxMappedInstanceID = 0;

        for (const auto& record : records)
        {
            maxMappedInstanceID =
                std::max(
                    maxMappedInstanceID,
                    record.instanceID
                );
        }

        // Determine the original triangle count for every mapped instance.
        std::vector<uint32_t> triangleCounts(
            size_t(maxMappedInstanceID) + 1,
            0u
        );

        const auto triangleInstanceIDs =
            pScene->getGeometryInstanceIDsByType(
                Scene::GeometryType::TriangleMesh
            );

        for (uint32_t instanceID : triangleInstanceIDs)
        {
            if (instanceID > maxMappedInstanceID)
                continue;

            const auto& instance =
                pScene->getGeometryInstance(instanceID);

            const MeshID meshID(instance.geometryID);
            const auto& mesh = pScene->getMesh(meshID);

            triangleCounts[instanceID] =
                mesh.indexCount / 3;
        }

        std::vector<uint2> instanceRanges(
            size_t(maxMappedInstanceID) + 1
        );

        uint64_t totalTriangleSlots64 = 0;

        for (uint32_t instanceID = 0;
            instanceID <= maxMappedInstanceID;
            ++instanceID)
        {
            if (totalTriangleSlots64 >
                uint64_t(std::numeric_limits<uint32_t>::max()))
            {
                FALCOR_THROW(
                    "Bistro runtime lightmap lookup exceeds uint32 indexing."
                );
            }

            instanceRanges[instanceID] =
                uint2(
                    uint32_t(totalTriangleSlots64),
                    triangleCounts[instanceID]
                );

            totalTriangleSlots64 +=
                triangleCounts[instanceID];
        }

        if (totalTriangleSlots64 >
            uint64_t(std::numeric_limits<uint32_t>::max()))
        {
            FALCOR_THROW(
                "Bistro runtime lightmap lookup exceeds uint32 indexing."
            );
        }

        const uint32_t totalTriangleSlots =
            uint32_t(totalTriangleSlots64);

        std::vector<uint32_t> pageIndices(
            totalTriangleSlots,
            0xffffffffu
        );

        std::vector<float4> uv01(
            totalTriangleSlots,
            float4(0.f)
        );

        std::vector<float2> uv2(
            totalTriangleSlots,
            float2(0.f)
        );

        uint64_t loadedRecordCount = 0;

        for (const auto& record : records)
        {
            if (record.instanceID > maxMappedInstanceID)
                continue;

            const uint2 range =
                instanceRanges[record.instanceID];

            if (record.triangleID >= range.y)
            {
                FALCOR_THROW(
                    "Bistro atlas mapping has invalid triangleID={} "
                    "for instanceID={} (triangleCount={}).",
                    record.triangleID,
                    record.instanceID,
                    range.y
                );
            }

            if (record.pageIndex >= header.pageCount)
            {
                FALCOR_THROW(
                    "Bistro atlas mapping has invalid pageIndex={}.",
                    record.pageIndex
                );
            }

            const uint32_t lookupIndex =
                range.x + record.triangleID;

            pageIndices[lookupIndex] =
                record.pageIndex;

            uv01[lookupIndex] =
                float4(
                    record.uv[0],
                    record.uv[1],
                    record.uv[2],
                    record.uv[3]
                );

            uv2[lookupIndex] =
                float2(
                    record.uv[4],
                    record.uv[5]
                );

            ++loadedRecordCount;
        }

        if (instanceRanges.empty() ||
            pageIndices.empty())
        {
            FALCOR_THROW(
                "Bistro atlas mapping produced an empty runtime lookup."
            );
        }

        spBistroInstanceRanges =
            pDevice->createStructuredBuffer(
                sizeof(uint2),
                uint32_t(instanceRanges.size()),
                ResourceBindFlags::ShaderResource
            );

        spBistroInstanceRanges->setBlob(
            instanceRanges.data(),
            0,
            instanceRanges.size() * sizeof(uint2)
        );

        spBistroPageIndices =
            pDevice->createStructuredBuffer(
                sizeof(uint32_t),
                totalTriangleSlots,
                ResourceBindFlags::ShaderResource
            );

        spBistroPageIndices->setBlob(
            pageIndices.data(),
            0,
            pageIndices.size() * sizeof(uint32_t)
        );

        spBistroUV01 =
            pDevice->createStructuredBuffer(
                sizeof(float4),
                totalTriangleSlots,
                ResourceBindFlags::ShaderResource
            );

        spBistroUV01->setBlob(
            uv01.data(),
            0,
            uv01.size() * sizeof(float4)
        );

        spBistroUV2 =
            pDevice->createStructuredBuffer(
                sizeof(float2),
                totalTriangleSlots,
                ResourceBindFlags::ShaderResource
            );

        spBistroUV2->setBlob(
            uv2.data(),
            0,
            uv2.size() * sizeof(float2)
        );

        const std::vector<std::string> pageFiles =
            readBistroAtlasPageFiles(header.pageCount);

        for (uint32_t pageIndex = 0; pageIndex < header.pageCount; ++pageIndex)
        {
            auto& page = spBistroAtlasPages[pageIndex];
            page = Texture::createFromFile(
                pDevice, pageFiles[pageIndex], false, false, ResourceBindFlags::ShaderResource);
            if (!page)
                FALCOR_THROW("Failed to load atlas page {} from '{}'.", pageIndex, pageFiles[pageIndex]);
            if (page->getWidth() != header.width || page->getHeight() != header.height)
                FALCOR_THROW("Atlas page {} dimensions do not match the mapping.", pageIndex);
            page->setName("BistroXAtlasPage" + std::to_string(pageIndex));
        }
        // Keep every descriptor valid when the scene uses fewer than ten pages.
        for (uint32_t pageIndex = header.pageCount; pageIndex < kBistroTestMaxPages; ++pageIndex)
            spBistroAtlasPages[pageIndex] = spBistroAtlasPages[0];

        sBistroPageCount =
            header.pageCount;

        sBistroMaxInstanceID =
            maxMappedInstanceID;

        sBistroLightmapEnabled = true;

        logInfo(
            "Loaded Bistro xatlas lightmap test: "
            "pages={} size={}x{} records={} "
            "maxInstanceID={} lookupSlots={}.",
            header.pageCount,
            header.width,
            header.height,
            loadedRecordCount,
            sBistroMaxInstanceID,
            totalTriangleSlots
        );
    }

    void bindBistroXAtlasLightmapTest(
        ShaderVar var,
        const ref<Sampler>& pLinearSampler
    )
    {
        var["gLinearSampler"] =
            pLinearSampler;

        var["BistroLightmapCB"]["gBistroLightmapEnabled"] =
            sBistroLightmapEnabled ? 1u : 0u;

        var["BistroLightmapCB"]["gBistroPageCount"] =
            sBistroPageCount;

        var["BistroLightmapCB"]["gBistroMaxInstanceID"] =
            sBistroMaxInstanceID;

        var["BistroLightmapCB"]["gBistroPadding"] =
            0u;

        if (!sBistroLightmapEnabled)
            return;

        for (uint32_t pageIndex = 0; pageIndex < kBistroTestMaxPages; ++pageIndex)
            var["gBistroAtlasPages"][pageIndex] = spBistroAtlasPages[pageIndex];

        var["gBistroInstanceRanges"] =
            spBistroInstanceRanges;

        var["gBistroPageIndices"] =
            spBistroPageIndices;

        var["gBistroUV01"] =
            spBistroUV01;

        var["gBistroUV2"] =
            spBistroUV2;
    }
}

void AdaptiveSHDemo::captureCurrentCameraOrbit()
{
    if (!mpScene)
    {
        return;
    }

    const auto pCamera = mpScene->getCamera();
    if (!pCamera)
    {
        return;
    }

    constexpr float kEpsilon = 1e-4f;

    if (
        mCameraOrbitRadiusX < kEpsilon ||
        mCameraOrbitRadiusZ < kEpsilon
        )
    {
        logWarning(
            "Cannot capture ellipse orbit with a zero radius."
        );
        return;
    }

    const float3 position =
        pCamera->getPosition();

    const float3 worldOffset =
        position - mCameraOrbitCenter;

    mCameraOrbitHeight = worldOffset.y;

    // Transform the camera offset into the ellipse's local frame.
    const float yaw =
        degreesToRadians(mCameraOrbitYawDeg);

    const float cosYaw = std::cos(yaw);
    const float sinYaw = std::sin(yaw);

    const float localX =
        cosYaw * worldOffset.x +
        sinYaw * worldOffset.z;

    const float localZ =
        -sinYaw * worldOffset.x +
        cosYaw * worldOffset.z;

    // Uniformly scale the ellipse so the current camera position lies on it,
    // while preserving the X/Z radius ratio.
    const float normalizedX =
        localX / mCameraOrbitRadiusX;

    const float normalizedZ =
        localZ / mCameraOrbitRadiusZ;

    const float ellipseScale =
        std::sqrt(
            normalizedX * normalizedX +
            normalizedZ * normalizedZ
        );

    if (ellipseScale < kEpsilon)
    {
        logWarning(
            "Cannot capture orbit: camera is too close "
            "to the orbit center."
        );
        return;
    }

    mCameraOrbitRadiusX *= ellipseScale;
    mCameraOrbitRadiusZ *= ellipseScale;

    mCameraOrbitAngle = std::atan2(
        localZ / mCameraOrbitRadiusZ,
        localX / mCameraOrbitRadiusX
    );
}

void AdaptiveSHDemo::startCameraOrbit()
{
    if (!mpScene || !mpScene->getCamera())
    {
        logWarning("Cannot start camera orbit: no camera is available.");
        return;
    }

    if (mCameraOrbitUseCurrentPoseOnStart)
    {
        captureCurrentCameraOrbit();
    }
    else
    {
        applyCameraOrbitPose();
    }

    mCameraOrbitLastTime = std::chrono::steady_clock::now();
    mCameraOrbitActive = true;
}

void AdaptiveSHDemo::stopCameraOrbit()
{
    mCameraOrbitActive = false;
}

void AdaptiveSHDemo::toggleCameraOrbit()
{
    if (mCameraOrbitActive)
    {
        stopCameraOrbit();
    }
    else
    {
        startCameraOrbit();
    }
}

void AdaptiveSHDemo::updateCameraOrbit()
{
    if (!mCameraOrbitActive || !mpScene)
    {
        return;
    }

    const auto currentTime = std::chrono::steady_clock::now();

    float deltaTime = std::chrono::duration<float>(
        currentTime - mCameraOrbitLastTime
    ).count();

    mCameraOrbitLastTime = currentTime;

    // Avoid jumps after debugging or a long frame.
    deltaTime = std::clamp(deltaTime, 0.0f, 0.1f);

    mCameraOrbitAngle +=
        degreesToRadians(mCameraOrbitSpeedDeg) * deltaTime;

    // Prevent the angle from growing indefinitely.
    if (mCameraOrbitAngle > 2.0f * kPi)
    {
        mCameraOrbitAngle -= 2.0f * kPi;
    }
    else if (mCameraOrbitAngle < -2.0f * kPi)
    {
        mCameraOrbitAngle += 2.0f * kPi;
    }

    applyCameraOrbitPose();
}

void AdaptiveSHDemo::applyCameraOrbitPose()
{
    if (!mpScene)
    {
        return;
    }

    const auto pCamera = mpScene->getCamera();
    if (!pCamera)
    {
        return;
    }

    const float theta = mCameraOrbitAngle;
    const float yaw =
        degreesToRadians(mCameraOrbitYawDeg);

    // Ellipse in its local XZ plane.
    const float localX =
        mCameraOrbitRadiusX * std::cos(theta);

    const float localZ =
        mCameraOrbitRadiusZ * std::sin(theta);

    // Rotate the ellipse around world Y.
    const float cosYaw = std::cos(yaw);
    const float sinYaw = std::sin(yaw);

    const float offsetX =
        cosYaw * localX -
        sinYaw * localZ;

    const float offsetZ =
        sinYaw * localX +
        cosYaw * localZ;

    const float3 cameraPosition(
        mCameraOrbitCenter.x + offsetX,
        mCameraOrbitCenter.y + mCameraOrbitHeight,
        mCameraOrbitCenter.z + offsetZ
    );

    //const float3 lookAtPoint =
    //    mCameraOrbitCenter + mCameraLookAtOffset;
    pCamera->setPosition(cameraPosition);
    pCamera->setTarget(mCameraLookAtTarget);

    //pCamera->setPosition(cameraPosition);
    //pCamera->setTarget(lookAtPoint);
}

Properties AdaptiveSHDemo::getProperties() const
{
    Properties props;
    props["runtimeGlass"] = mEnableGlass;
    props["shInterpolation"] = mShInterpolation;
    props["shSeam"] = mShSeam;
    props["shNonNegative"] = mShNonNegative;
    props["glassRefraction"] = mGlassRefractionStrength;
    props["liquidOpticalDepth"] = mLiquidOpticalDepth;
    return props;
}

RenderPassReflection AdaptiveSHDemo::reflect(const CompileData& compileData)
{
    // Define the required resources here
    RenderPassReflection reflector;
    const uint2 sz = RenderPassHelpers::calculateIOSize(mOutputSizeSelection, mFixedOutputSize, compileData.defaultTexDims);
    // REMARK MSAA is set via texture sample count. Note that all fbo attachment have to have same sample count.
    reflector.addOutput("output", "Color").texture2D(sz.x, sz.y, 4).format(ResourceFormat::RGBA32Float)
        .bindFlags(ResourceBindFlags::RenderTarget | ResourceBindFlags::ShaderResource);
    //reflector.addOutput("output", "Color").texture2D(mLightmapWidth, mLightmapHeight, 4).format(ResourceFormat::RGBA32Float);
    reflector.addOutput("depth", "Depth buffer")
        .format(ResourceFormat::D32Float)
        .bindFlags(ResourceBindFlags::DepthStencil | ResourceBindFlags::ShaderResource)
        .texture2D(sz.x, sz.y, 4);

    reflector.addInternal("dynamicRaw", "Dynamic raw lighting")
        .format(ResourceFormat::RGBA32Float)
        .bindFlags(ResourceBindFlags::RenderTarget | ResourceBindFlags::ShaderResource)
        .texture2D(sz.x, sz.y);

    reflector.addOutput("dynamicFiltered", "Dynamic filtered lighting")
        .format(ResourceFormat::RGBA32Float)
        .bindFlags(ResourceBindFlags::UnorderedAccess | ResourceBindFlags::ShaderResource)
        .texture2D(sz.x, sz.y);

    reflector.addInternal("dynamicMask", "Dynamic mask")
        .format(ResourceFormat::R8Unorm)
        .bindFlags(ResourceBindFlags::RenderTarget | ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess)
        .texture2D(sz.x, sz.y);

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
        updateCameraOrbit();

        // Runtime measurement: the timers of this ring slot were written kTimerRing frames ago; read them before reuse.
        if (mFrameTimers.empty())
        {
            mFrameTimers.resize(kTimerRing);
            for (FrameTimers& t : mFrameTimers)
            {
                t.total = GpuTimer::create(mpDevice);
                t.shading = GpuTimer::create(mpDevice);
                t.glass = GpuTimer::create(mpDevice);
            }
        }
        FrameTimers& timers = mFrameTimers[mTimerFrame % kTimerRing];
        if (timers.measured)
        {
            mMeasGpuTotalMs.push_back(timers.total->getElapsedTime());
            mMeasGpuShadingMs.push_back(timers.shading->getElapsedTime());
            mMeasGpuGlassMs.push_back(timers.glassUsed ? timers.glass->getElapsedTime() : 0.0);
            timers.measured = false;
            if (!mMeasureActive && mMeasGpuTotalMs.size() == mMeasCpuMs.size()) writeRuntimeMeasurement();
        }
        if (mMeasureRequested && !mMeasureActive && mMeasGpuTotalMs.size() == mMeasCpuMs.size())
        {
            mMeasureRequested = false;
            mMeasureActive = true;
            mMeasureFrame = 0;
            mMeasCpuMs.clear(); mMeasGpuTotalMs.clear(); mMeasGpuShadingMs.clear(); mMeasGpuGlassMs.clear();
            mLastMeasurement = "running...";
        }
        const bool measureThisFrame = mMeasureActive && mMeasureFrame >= kMeasureWarmupFrames;
        timers.total->begin();
        timers.shading->begin();

        // ------------------------------------------------------------------
        // PASS 1: STATIC GEOMETRY (lightmaps)
        // ------------------------------------------------------------------
        bindShadingVars(mpStaticVars->getRootVar());
        mpScene->rasterize(pRenderContext, mpGraphicsState.get(), mpStaticVars.get(), mpRasterState, mpRasterState);
        timers.shading->end();
        timers.glassUsed = mEnableGlass;
        if (mEnableGlass) timers.glass->begin();

        if (mEnableGlass)
        {
            const uint32_t width = pTargetFbo->getWidth(), height = pTargetFbo->getHeight();
            if (!mpOpaqueColor || mpOpaqueColor->getWidth() != width || mpOpaqueColor->getHeight() != height)
            {
                auto flags = ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess;
                mpOpaqueColor = mpDevice->createTexture2D(width, height, ResourceFormat::RGBA32Float, 1, 1, nullptr, flags);
                mpOpaqueDepth = mpDevice->createTexture2D(width, height, ResourceFormat::R32Float, 1, 1, nullptr, flags);
            }
            for (uint32_t stage = mHasLiquidMaterials ? 0u : 1u; stage < 2; ++stage)
            {
                auto resolve = mpOpaqueResolvePass->getRootVar();
                resolve["gColorMS"] = pTargetFbo;
                resolve["gDepthMS"] = pDepth;
                resolve["gColor"] = mpOpaqueColor;
                resolve["gDepth"] = mpOpaqueDepth;
                mpOpaqueResolvePass->execute(pRenderContext, width, height);
                // Release attachment SRVs before they are rebound for glass drawing.
                resolve["gColorMS"] = ref<Texture>(); resolve["gDepthMS"] = ref<Texture>();
                auto glass = mpGlassVars->getRootVar();
                glass["gLiquidMaterials"] = mpLiquidMaterialMask;
                glass["GlassCB"]["gLiquidStage"] = stage == 0 ? 1u : 0u;
                glass["GlassCB"]["gLiquidOpticalDepth"] = mLiquidOpticalDepth;
                bindBistroXAtlasLightmapTest(glass, mpLinearSampler);
                glass["gOpaqueColor"] = mpOpaqueColor;
                glass["gOpaqueDepth"] = mpOpaqueDepth;
                glass["GlassCB"]["gHasEnvironment"] = mpScene->useEnvLight() ? 1u : 0u;
                glass["GlassCB"]["gSceneScale"] = std::max(mpScene->getSceneBounds().radius(), 1e-3f);
                glass["GlassCB"]["gRefractionStrength"] = mGlassRefractionStrength;
                // Stage 0 composites liquid against opaque surfaces. Stage 1
                // resolves that result again so outer glass sees the liquid.
                // Depth writes keep the nearest surface within each stage.
                mpScene->rasterize(pRenderContext, mpGlassState.get(), mpGlassVars.get(), mpRasterState, mpRasterState);
            }
        }

        if (mEnableGlass) timers.glass->end();
        timers.total->end();
        timers.total->resolve();
        timers.shading->resolve();
        if (mEnableGlass) timers.glass->resolve();
        timers.measured = measureThisFrame;
        mTimerFrame++;

        // CPU frame time: from this point in the previous frame to this point now (the whole Mogwai frame).
        const auto now = std::chrono::high_resolution_clock::now();
        const double frameMs = mHasLastFrameTime ? std::chrono::duration<double, std::milli>(now - mLastFrameTime).count() : 0.0;
        mLastFrameTime = now;
        mHasLastFrameTime = true;
        if (mMeasureActive)
        {
            if (measureThisFrame) mMeasCpuMs.push_back(frameMs);
            if (++mMeasureFrame >= kMeasureWarmupFrames + kMeasureFrames) mMeasureActive = false; // GPU samples follow
        }

        // View dump (UI): after the frame, before the grid overlays, so the saved frame shows the scene only.
        if (mDumpViewRequested)
        {
            mDumpViewRequested = false;
            dumpView(pRenderContext, pTargetFbo);
        }

        if (mbShowAdaptiveGrid)
        {
            mpProbeVisualizePass->setCameraData(
                mpScene->getCamera()->getViewProjMatrix()
            );

            mpProbeVisualizePass->execute(pRenderContext, mpFbo);
        }
        if (mbShowFDPlane && mpFDPlaneVisPass)
        {
            mpFDPlaneVisPass->setCameraData(mpScene->getCamera()->getViewProjMatrix());
            mpFDPlaneVisPass->execute(pRenderContext, mpFbo);
        }
    }
}
namespace
{
void meanStd(const std::vector<double>& v, double& mean, double& sd)
{
    mean = 0.0; sd = 0.0;
    if (v.empty()) return;
    for (double x : v) mean += x;
    mean /= double(v.size());
    for (double x : v) sd += (x - mean) * (x - mean);
    sd = std::sqrt(sd / double(v.size()));
}
} // namespace

// One row per measurement in AdaptiveSHDemo_RuntimeTiming.csv (new file; the old AdaptiveSHDemo_RuntimeFPS.csv is left as is).
void AdaptiveSHDemo::writeRuntimeMeasurement()
{
    double cpuMean, cpuSd, totMean, totSd, shMean, shSd, glMean, glSd;
    meanStd(mMeasCpuMs, cpuMean, cpuSd);
    meanStd(mMeasGpuTotalMs, totMean, totSd);
    meanStd(mMeasGpuShadingMs, shMean, shSd);
    meanStd(mMeasGpuGlassMs, glMean, glSd);

    uint64_t probeCount = 0;
#if CURRENT_PROBE_MODE == PROBE_MODE_ADAPTIVE
    const std::string modeName = "Adaptive";
    if (mAdaptiveProbeVolume) probeCount = mAdaptiveProbeVolume->getProbes().size();
    const std::string seamName = mShSeam == 1 ? "Oat" : (mShSeam == 2 ? "Constrained" : "None");
#else
    const std::string modeName = "Uniform";
    if (mUniformProbeVolume)
    {
        const uint3 res = mUniformProbeVolume->getCellResolution();
        probeCount = uint64_t(res.x + 1) * uint64_t(res.y + 1) * uint64_t(res.z + 1);
    }
    const std::string seamName = "-";
#endif
    const std::string interpName = mShInterpolation == 1 ? "Trilinear" : (mShInterpolation == 2 ? "LimitedHermite" : "Hermite");
    std::string sceneName = "-";
    if (mpScene)
    {
        const auto paths = mpScene->getImportPaths();
        if (!paths.empty()) sceneName = paths.front().stem().string();
    }
    uint32_t width = 0, height = 0;
    if (mpFbo && mpFbo->getColorTexture(0))
    {
        width = mpFbo->getColorTexture(0)->getWidth();
        height = mpFbo->getColorTexture(0)->getHeight();
    }

    const char* kFile = "AdaptiveSHDemo_RuntimeTiming.csv";
    const bool newFile = !std::ifstream(kFile).good();
    std::ofstream out(kFile, std::ios::app);
    if (newFile)
        out << "Scene,Mode,ProbeFile,ProbeCount,Width,Height,Interpolation,Seam,Deringing,Glass,WarmupFrames,MeasuredFrames,"
               "FrameMs,FrameStdMs,FPS,GpuTotalMs,GpuTotalStdMs,GpuShadingMs,GpuShadingStdMs,GpuGlassMs,GpuGlassStdMs\n";
    out << sceneName << "," << modeName << "," << loadFromFileName << "," << probeCount << "," << width << "," << height << ","
        << interpName << "," << seamName << "," << (mShNonNegative ? 1 : 0) << "," << (mEnableGlass ? 1 : 0) << ","
        << kMeasureWarmupFrames << "," << mMeasCpuMs.size() << "," << std::fixed << std::setprecision(4)
        << cpuMean << "," << cpuSd << "," << (cpuMean > 0.0 ? 1000.0 / cpuMean : 0.0) << ","
        << totMean << "," << totSd << "," << shMean << "," << shSd << "," << glMean << "," << glSd << "\n";

    char buf[256];
    std::snprintf(buf, sizeof(buf), "frame %.3f ms (%.1f fps), GPU total %.3f ms, shading %.3f ms, glass %.3f ms",
                  cpuMean, cpuMean > 0.0 ? 1000.0 / cpuMean : 0.0, totMean, shMean, glMean);
    mLastMeasurement = buf;
    logInfo("Runtime measurement ({}, {}, deringing {}): {}", sceneName, loadFromFileName, mShNonNegative ? "on" : "off", buf);
}

bool AdaptiveSHDemo::onKeyEvent(
    const KeyboardEvent& keyEvent
)
{
    if (
        keyEvent.type == KeyboardEvent::Type::KeyPressed &&
        keyEvent.key == Input::Key::O
        )
    {
        toggleCameraOrbit();
        return true;
    }

    return false;
}

// The probes of every enabled FD area (PrecomputeSHCoefficients/FDPlaneSettings.h) as small cubes in
// the area's colour, 40% of the probe spacing wide (8 mm for a 2 m area at 100 x 100).
void AdaptiveSHDemo::updateFDPlaneVis()
{
    if (!mpScene || !mpFDPlaneVisPass) return;
    const uint32_t res = mbFDPlaneRetraceRes ? fdplane::kRetraceRes : fdplane::kRes;
    std::vector<ProbeVisualizePass::PointGroup> groups;
    for (const fdplane::Area& a : fdplane::kAreas)
    {
        if (!a.enabled) continue;
        groups.push_back({ fdplane::positions(a, res), 0.2f * fdplane::spacing(a, res), a.color });
        int u, v;
        fdplane::planeAxes(a.normal, u, v);
        const float n = a.centre[3 - u - v];
        logInfo("FD area {}: {} x {} probes at {} = {}, {} in [{}, {}], {} in [{}, {}].", a.name, res, res,
                fdplane::normalName(a.normal), n, fdplane::axisName(u), a.centre[u] - a.halfSize, a.centre[u] + a.halfSize,
                fdplane::axisName(v), a.centre[v] - a.halfSize, a.centre[v] + a.halfSize);
    }
    mpFDPlaneVisPass->setPointGroups(groups);
}

// Lightmap, glass flag and probe grid of the static shading program (psMain); the dump program uses the same.
void AdaptiveSHDemo::bindShadingVars(ShaderVar applyVar)
{
    applyVar["GlassCB"]["gGlassEnabled"] = mEnableGlass ? 1u : 0u;
    bindBistroXAtlasLightmapTest(
        applyVar,
        mpLinearSampler
    );
    // Atlas-mapped geometry uses lightmaps; unmapped geometry uses the selected SH grid.
#if CURRENT_PROBE_MODE == PROBE_MODE_ADAPTIVE
    applyVar["gCornerBuffer"] = mAdaptiveProbeVolume->getCornerBuffer();
    applyVar["gProbeBuffer"] = mAdaptiveProbeVolume->getProbeBuffer();

    applyVar["gSeedProbeIndices"] = mAdaptiveProbeVolume->getSeedProbeIndexBuffer();

    // The uploaded seed (a grid without a seed gets one 1 x 1 x 1 cell): the shader derives node bounds from it.
    applyVar["SeedGridCB"]["gUseSeedGrid"] = 1u;
    applyVar["SeedGridCB"]["gSeedMinPoint"] = mAdaptiveProbeVolume->getGpuSeedMinPoint();
    applyVar["SeedGridCB"]["gSeedCellSize"] = mAdaptiveProbeVolume->getGpuSeedCellSize();
    applyVar["SeedGridCB"]["gSeedResolution"] = mAdaptiveProbeVolume->getGpuSeedResolution();
    applyVar["SeedGridCB"]["gLatticeCellSize"] = mAdaptiveProbeVolume->getGpuLatticeCellSize();
#else
    mUniformProbeVolume->bindShaderData(applyVar);
#endif
}

// Writes <scene tag>_View_<name>.dump (binary, little-endian) in the working directory:
//   char[8] "SHVIEW1", uint32 width, uint32 height,
//   float[16] camera: position xyz, target xyz, up xyz, focal length, frame height, aspect ratio, near, far, 0, 0,
//   then 5 planes of width x height float4, row 0 = top: position (w = 1 where the SH grid shades the pixel),
//   shading normal, albedo, emission, and the demo frame (MSAA resolved; what this pass rendered this frame).
// Plus <...>.camera.txt with the camera as Mogwai Python, to restore the view. Pixels are single-sample.
void AdaptiveSHDemo::dumpView(RenderContext* pRenderContext, const ref<Texture>& pOutput)
{
    if (!mpScene || !mpDumpProgram) return;
    const uint32_t w = pOutput->getWidth(), h = pOutput->getHeight();
    if (!mpDumpFbo || mpDumpFbo->getWidth() != w || mpDumpFbo->getHeight() != h)
    {
        mpDumpFbo = Fbo::create(mpDevice);
        for (uint32_t i = 0; i < 4; ++i)
            mpDumpFbo->attachColorTarget(mpDevice->createTexture2D(w, h, ResourceFormat::RGBA32Float, 1, 1, nullptr,
                ResourceBindFlags::RenderTarget | ResourceBindFlags::ShaderResource), i);
        mpDumpFbo->attachDepthStencilTarget(mpDevice->createTexture2D(w, h, ResourceFormat::D32Float, 1, 1, nullptr,
            ResourceBindFlags::DepthStencil));
    }
    pRenderContext->clearFbo(mpDumpFbo.get(), float4(0.f), 1.f, 0, FboAttachmentType::All);
    mpDumpState->setFbo(mpDumpFbo);
    bindShadingVars(mpDumpVars->getRootVar());
    mpScene->rasterize(pRenderContext, mpDumpState.get(), mpDumpVars.get(), mpRasterState, mpRasterState);

    ref<Texture> pFrame = mpDevice->createTexture2D(w, h, ResourceFormat::RGBA32Float, 1, 1, nullptr,
        ResourceBindFlags::RenderTarget | ResourceBindFlags::ShaderResource);
    if (pOutput->getSampleCount() > 1) pRenderContext->resolveResource(pOutput, pFrame);
    else pRenderContext->copyResource(pFrame.get(), pOutput.get());

    std::vector<std::vector<uint8_t>> planes;
    for (uint32_t i = 0; i < 4; ++i) planes.push_back(pRenderContext->readTextureSubresource(mpDumpFbo->getColorTexture(i).get(), 0));
    planes.push_back(pRenderContext->readTextureSubresource(pFrame.get(), 0));

    const std::string stem = (mViewDumpScene.empty() ? std::string("Scene") : mViewDumpScene) + "_View_" +
                             (mViewDumpName.empty() ? std::string("view") : mViewDumpName);
    const std::string file = stem + ".dump";
    std::ofstream out(file, std::ios::binary);
    if (!out) { logError("[View dump] Cannot write {}.", file); return; }
    const char magic[8] = { 'S', 'H', 'V', 'I', 'E', 'W', '1', 0 };
    out.write(magic, 8);
    out.write(reinterpret_cast<const char*>(&w), 4);
    out.write(reinterpret_cast<const char*>(&h), 4);
    const auto& cam = mpScene->getCamera();
    const float3 p = cam->getPosition(), t = cam->getTarget(), u = cam->getUpVector();
    const float camData[16] = { p.x, p.y, p.z, t.x, t.y, t.z, u.x, u.y, u.z, cam->getFocalLength(), cam->getFrameHeight(),
                                cam->getAspectRatio(), cam->getNearPlane(), cam->getFarPlane(), 0.f, 0.f };
    out.write(reinterpret_cast<const char*>(camData), sizeof(camData));
    const size_t planeBytes = size_t(w) * h * 16;
    for (const auto& pl : planes)
    {
        if (pl.size() < planeBytes) { logError("[View dump] Read-back of {} bytes, expected {}.", pl.size(), planeBytes); return; }
        out.write(reinterpret_cast<const char*>(pl.data()), std::streamsize(planeBytes));
    }
    out.close();
    size_t shPixels = 0;
    const float* posPlane = reinterpret_cast<const float*>(planes[0].data());
    for (size_t i = 0; i < size_t(w) * h; ++i) shPixels += posPlane[4 * i + 3] > 0.5f;
    std::ofstream camOut(stem + ".camera.txt");
    camOut << std::setprecision(9) << "# Mogwai: paste into the script console to restore the view of " << file << "\n"
           << "camera = m.scene.camera\n"
           << "camera.position = float3(" << p.x << ", " << p.y << ", " << p.z << ")\n"
           << "camera.target = float3(" << t.x << ", " << t.y << ", " << t.z << ")\n"
           << "camera.up = float3(" << u.x << ", " << u.y << ", " << u.z << ")\n"
           << "camera.focalLength = " << cam->getFocalLength() << "\n"
           << "# frame " << w << " x " << h << ", grid file in the demo: " << loadFromFileName << "\n";
    mLastViewDump = file + " (" + std::to_string(w) + " x " + std::to_string(h) + ", " + std::to_string(shPixels) + " SH pixels)";
    logInfo("[View dump] {}: {} x {}, {} pixels shaded by the SH grid.", file, w, h, shPixels);
}

void AdaptiveSHDemo::loadAdaptiveGrid()
{
#if CURRENT_PROBE_MODE == PROBE_MODE_ADAPTIVE
    mAdaptiveProbeVolume = AdaptiveProbeVolume::create(mpDevice);

    // A missing grid file leaves the volume empty, and the upload then fails with "Can't create GPU buffer of size
    // zero". Stop here with the file name instead.
    if (!std::ifstream(loadFromFileName))
        FALCOR_THROW("AdaptiveSHDemo: grid file '{}' not found in the working directory (loadFromFileName).", loadFromFileName);
    mAdaptiveProbeVolume->loadFromFile(loadFromFileName);
    mAdaptiveProbeVolume->setNonNegativeIrradiance(mShNonNegative);
    mAdaptiveProbeVolume->setHangingConstraint(mShSeam == 2 ? int(mShInterpolation) : -1);
    mAdaptiveProbeVolume->uploadToGPU();
    if (mpProbeVisualizePass) mpProbeVisualizePass->setVolumeData(mAdaptiveProbeVolume->getProbes());
#endif
}

void AdaptiveSHDemo::loadUniformGrid()
{
#if CURRENT_PROBE_MODE != PROBE_MODE_ADAPTIVE
    mUniformProbeVolume = UniformProbeVolume::create(mpDevice);
    if (!std::ifstream(loadFromFileName))
        FALCOR_THROW("AdaptiveSHDemo: grid file '{}' not found in the working directory (loadFromFileName).", loadFromFileName);
    mUniformProbeVolume->loadFromFile(loadFromFileName);
    mUniformProbeVolume->setNonNegativeIrradiance(mShNonNegative);
    mUniformProbeVolume->uploadToGPU();
    if (mpProbeVisualizePass)
        mpProbeVisualizePass->setUniformVolumeData(
            mUniformProbeVolume->getMinPoint(), mUniformProbeVolume->getCellSize(), mUniformProbeVolume->getCellResolution());
#endif
}

DefineList AdaptiveSHDemo::getShInterpolationDefines() const
{
    DefineList d;
    d.add("SH_INTERP", std::to_string(mShInterpolation));
    d.add("SH_SEAM_OAT", mShSeam == 1 ? "1" : "0");
    d.add("SH_NONNEG", mShNonNegative ? "1" : "0");
    return d;
}

// Changing the defines recompiles the programs; the vars are rebuilt for the new reflection (all of them are bound
// every frame in execute()).
void AdaptiveSHDemo::applyShInterpolationDefines()
{
    const DefineList d = getShInterpolationDefines();
    if (mpStaticProgram)
    {
        mpStaticProgram->addDefines(d);
        mpStaticVars = ProgramVars::create(mpDevice, mpStaticProgram->getReflector());
    }
    if (mpGlassProgram)
    {
        mpGlassProgram->addDefines(d);
        mpGlassVars = ProgramVars::create(mpDevice, mpGlassProgram->getReflector());
    }
}

void AdaptiveSHDemo::renderUI(Gui::Widgets& widget) {
    widget.checkbox("Runtime glass", mEnableGlass);
    widget.var("Glass refraction", mGlassRefractionStrength, 0.f, 3.f);
    widget.var("Liquid color depth", mLiquidOpticalDepth, 0.f, 10.0f, 0.001f);
    widget.text("Loaded probe file: " + loadFromFileName);
    {
        // Same interpolation modes as the placement evaluation (PlacementComparison.cpp).
        static const Gui::DropdownList kShInterpolationList = {
            { 0, "Hermite" }, { 1, "Trilinear" }, { 2, "Hermite, limited slopes" } };
        const uint32_t oldInterp = mShInterpolation, oldSeam = mShSeam;
        bool changed = widget.dropdown("SH interpolation", kShInterpolationList, mShInterpolation);
#if CURRENT_PROBE_MODE == PROBE_MODE_ADAPTIVE
        static const Gui::DropdownList kSeamList = {
            { 0, "None" }, { 1, "Oat (runtime lookups)" }, { 2, "Hanging-node constraint (baked)" } };
        changed |= widget.dropdown("Seam handling", kSeamList, mShSeam);
        widget.tooltip("Oat: lattice interpolation next to finer leaves, in the shader (slow). Constraint: hanging corners take "
                       "the coarse side's values when the grid is loaded; the shader runs the plain interpolation.");
#endif
        if (changed) applyShInterpolationDefines();
#if CURRENT_PROBE_MODE == PROBE_MODE_ADAPTIVE
        // The constraint is baked into the corners for one interpolation: load the grid again when it is switched on or
        // off, or when the interpolation changes while it is on.
        if ((mShSeam == 2) != (oldSeam == 2) || (mShSeam == 2 && mShInterpolation != oldInterp)) loadAdaptiveGrid();
#endif
        // Sloan deringing of the loaded probes (load the grid again) and of the interpolated SH at shading points, where
        // Hermite can overshoot between windowed probes (shader define SH_NONNEG).
        if (widget.checkbox("SH deringing (Sloan 2017)", mShNonNegative))
        {
            applyShInterpolationDefines();
#if CURRENT_PROBE_MODE == PROBE_MODE_ADAPTIVE
            loadAdaptiveGrid();
#else
            loadUniformGrid();
#endif
        }
        widget.tooltip("Order-2 SH irradiance goes below zero on normals facing away from a small, bright emitter close to "
                       "the probe; the clamp to zero then removes only the bright channels (coloured patches). On: Sloan's "
                       "deringing (SIGGRAPH Asia 2017): every such probe gets the least window that makes its irradiance "
                       "non-negative, the same for R, G and B (when the grid is loaded), and so does the interpolated SH at "
                       "each shading point.");
#if CURRENT_PROBE_MODE == PROBE_MODE_ADAPTIVE
        if (mAdaptiveProbeVolume && mShNonNegative)
            widget.text("Windowed corners: " + std::to_string(mAdaptiveProbeVolume->getNonNegativeCornerCount()));
#else
        if (mUniformProbeVolume && mShNonNegative)
            widget.text("Windowed probes: " + std::to_string(mUniformProbeVolume->getNonNegativeProbeCount()));
#endif
    }
    if (auto g = widget.group("Runtime measurement", true))
    {
        const std::string measureLabel = "Measure (" + std::to_string(kMeasureFrames) + " frames)";
        if (g.button(measureLabel.c_str()) && !mMeasureActive) mMeasureRequested = true;
        g.tooltip("Times the next frames after a short warm-up and appends a row to AdaptiveSHDemo_RuntimeTiming.csv with the "
                  "current settings. GPU times are this pass's own work and do not depend on v-sync; the frame time is the "
                  "whole frame, so turn v-sync off for it. Keep the camera still.");
        g.text("Last: " + (mLastMeasurement.empty() ? std::string("-") : mLastMeasurement));
    }
    if (auto g = widget.group("View dump (offline image evaluation)", false))
    {
        // The scene tag must be the precompute scene's .pyscene name (the evaluation reads <tag>_View_*.dump).
        g.textbox("Scene tag", mViewDumpScene);
        g.textbox("View name", mViewDumpName);
        if (g.button("Dump view")) mDumpViewRequested = true;
        if (!mLastViewDump.empty()) g.text("Last: " + mLastViewDump);
    }

    if (auto orbitGroup = widget.group("Camera Orbit", true))
    {
        orbitGroup.text("Hotkey: O");

        orbitGroup.text(
            mCameraOrbitActive
            ? "Status: running"
            : "Status: stopped"
        );
        // --------------------------------------------------------------
        // Position and target
        // --------------------------------------------------------------

        bool poseChanged = false;

        poseChanged |= orbitGroup.var(
            "Orbit center",
            mCameraOrbitCenter,
            -10.0f,
            10.0f,
            0.01f
        );

        //poseChanged |= orbitGroup.var(
        //    "Look-at offset",
        //    mCameraLookAtOffset,
        //    -5.0f,
        //    5.0f,
        //    0.01f
        //);

        poseChanged |= orbitGroup.var(
            "Look-at target",
            mCameraLookAtTarget,
            -10.0f,
            10.0f,
            0.01f
        );

        poseChanged |= orbitGroup.var(
            "Radius X",
            mCameraOrbitRadiusX,
            0.05f,
            20.0f,
            0.01f
        );

        poseChanged |= orbitGroup.var(
            "Radius Z",
            mCameraOrbitRadiusZ,
            0.05f,
            20.0f,
            0.01f
        );

        poseChanged |= orbitGroup.var(
            "Ellipse rotation (deg)",
            mCameraOrbitYawDeg,
            -180.0f,
            180.0f,
            0.5f
        );

        poseChanged |= orbitGroup.var(
            "Height",
            mCameraOrbitHeight,
            -10.0f,
            10.0f,
            0.01f
        );

        // Present the internal radian angle as degrees in the UI.
        float angleDeg = radiansToDegrees(mCameraOrbitAngle);

        if (
            orbitGroup.var(
                "Angle (deg)",
                angleDeg,
                -360.0f,
                360.0f,
                0.1f
            )
            )
        {
            mCameraOrbitAngle = degreesToRadians(angleDeg);
            poseChanged = true;
        }

        // Apply edits immediately, even while the animation is stopped.
        if (poseChanged)
        {
            applyCameraOrbitPose();
        }

        // --------------------------------------------------------------
        // Animation controls
        // --------------------------------------------------------------

        orbitGroup.var(
            "Speed (deg/s)",
            mCameraOrbitSpeedDeg,
            -90.0f,
            90.0f,
            0.1f
        );

        orbitGroup.checkbox(
            "Use current camera pose on start",
            mCameraOrbitUseCurrentPoseOnStart
        );

        if (orbitGroup.button("Capture current camera pose"))
        {
            captureCurrentCameraOrbit();
        }

        if (orbitGroup.button("Reverse direction"))
        {
            mCameraOrbitSpeedDeg = -mCameraOrbitSpeedDeg;
        }

        if (orbitGroup.button("Reset angle"))
        {
            mCameraOrbitAngle = 0.0f;
            applyCameraOrbitPose();
        }

        if (
            orbitGroup.button(
                mCameraOrbitActive
                ? "Stop orbit"
                : "Start orbit"
            )
            )
        {
            toggleCameraOrbit();
        }
    }

    if (widget.checkbox("Show SH Grid", mbShowAdaptiveGrid))
        requestRecompile();
    if (widget.checkbox("Show FD areas", mbShowFDPlane) && mbShowFDPlane) updateFDPlaneVis();
    const std::string retraceLabel = "FD areas: re-traced grid (" + std::to_string(fdplane::kRetraceRes) + " x " + std::to_string(fdplane::kRetraceRes) + ")";
    if (mbShowFDPlane && widget.checkbox(retraceLabel.c_str(), mbFDPlaneRetraceRes)) updateFDPlaneVis();
    // Level Visibility Controls
    if (mbShowAdaptiveGrid)
    {
        if (widget.checkbox(
            "Show Normal Voxels",
            mShowNormalVoxels
        ))
        {
            if (mpProbeVisualizePass)
            {
                mpProbeVisualizePass->setShowNormal(
                    mShowNormalVoxels
                );
            }
        }

        if (widget.checkbox(
            "Show Edge-Added Voxels",
            mbShowEdgeAddedVoxels
        ))
        {
            if (mpProbeVisualizePass)
            {
                mpProbeVisualizePass->setShowEdgeAdded(
                    mbShowEdgeAddedVoxels
                );
            }
        }

        if (widget.checkbox(
            "Draw Leaf Only",
            mbDrawLeafOnly
        ))
        {
            if (mpProbeVisualizePass)
            {
                mpProbeVisualizePass->setDrawLeafOnly(
                    mbDrawLeafOnly
                );
            }
        }

        if (auto g =
            widget.group("Octree Levels", true))
        {
            for (int i = 0; i < 8; ++i)
            {
                std::string label =
                    "Level " + std::to_string(i);

                if (g.checkbox(
                    label.c_str(),
                    mVisLevels[i]
                ))
                {
                    if (mpProbeVisualizePass)
                    {
                        mpProbeVisualizePass->toggleLevel(
                            i,
                            mVisLevels[i]
                        );
                    }
                }
            }
        }
    }
}

void AdaptiveSHDemo::setScene(RenderContext* pRenderContext, const ref<Scene>& pScene)
{
    mpScene = pScene;
    mMeasureRequested = true; // time the new scene once it has warmed up
    mpLiquidMaterialMask = nullptr;
    mHasLiquidMaterials = false;
    if (mpScene)
    {
        // Explicit roles from BistroInterior_Wine.pyscene. Keep containers such
        // as TransparentGlassWine in the outer-glass stage, not the liquid stage.
        // Edit this one list when another scene uses different liquid names.
        const std::array<std::string, 4> liquidNames = {"White_Wine", "Red_Wine", "Beer", "Water"};
        const auto& materials = mpScene->getMaterials();
        std::vector<uint32_t> liquidMask(std::max(size_t(1), materials.size()), 0u);
        for (size_t i = 0; i < materials.size(); ++i)
        {
            const auto& name = materials[i]->getName();
            if (std::find(liquidNames.begin(), liquidNames.end(), name) != liquidNames.end())
            {
                liquidMask[i] = 1u;
                mHasLiquidMaterials = true;
                logInfo("Runtime liquid material: '{}' (id={}).", name, i);
            }
        }
        mpLiquidMaterialMask = mpDevice->createStructuredBuffer(sizeof(uint32_t), uint32_t(liquidMask.size()), ResourceBindFlags::ShaderResource);
        mpLiquidMaterialMask->setBlob(liquidMask.data(), 0, liquidMask.size() * sizeof(uint32_t));
        std::ifstream check("AdaptiveSHDemo_RuntimeFPS.csv");
        if (!check.good())
        {
            std::ofstream out("AdaptiveSHDemo_RuntimeFPS.csv");
            out << "Mode,ProbeFile,ProbeCount,WarmupFrames,MeasuredFrames,MeanFrameMs,StdFrameMs,MeanFPS\n";
        }

        //init probe visual pass
        mpProbeVisualizePass = ProbeVisualizePass::create(mpDevice, mpScene->getSceneDefines());
        mpFDPlaneVisPass = ProbeVisualizePass::create(mpDevice, mpScene->getSceneDefines());
        updateFDPlaneVis();
        //Re-apply the visibility masks from our member variables
        for (int i = 0; i < 8; ++i)
        {
            mpProbeVisualizePass->toggleLevel(i, mVisLevels[i]);
        }
        // Restore Leaf Only
        mpProbeVisualizePass->setDrawLeafOnly(mbDrawLeafOnly);
        mpProbeVisualizePass->setShowEdgeAdded(
            mbShowEdgeAddedVoxels
        );
        mpProbeVisualizePass->setShowNormal(
            mShowNormalVoxels
        );
#if CURRENT_PROBE_MODE == PROBE_MODE_ADAPTIVE
        loadAdaptiveGrid();
#else
        loadUniformGrid();
#endif


        ProgramDesc previewDesc;
        previewDesc.addShaderModules(mpScene->getShaderModules());
        previewDesc.addShaderLibrary(kShaderFile)
            .vsEntry("vsMain")
            .psEntry("psMain");

        DefineList staticDefines = mpScene->getSceneDefines();
        staticDefines.add(getShInterpolationDefines());
        mpStaticProgram = Program::create(mpDevice, previewDesc, staticDefines);
        mpStaticVars = ProgramVars::create(mpDevice, mpStaticProgram->getReflector());

        // rasterizer state
        RasterizerState::Desc rasterDesc;
        rasterDesc.setFillMode(RasterizerState::FillMode::Solid);
        rasterDesc.setCullMode(RasterizerState::CullMode::None);
        //rasterDesc.setDepthBias(100000, 1.0f);
        mpRasterState = RasterizerState::create(rasterDesc);

        // default depth stencil state
        DepthStencilState::Desc dsDesc;
        ref<DepthStencilState> pDsState = DepthStencilState::create(dsDesc);

        mpGraphicsState = GraphicsState::create(mpDevice); 
        mpGraphicsState->setProgram(mpStaticProgram);
        mpGraphicsState->setRasterizerState(mpRasterState);
        mpGraphicsState->setFbo(mpFbo);
        mpGraphicsState->setDepthStencilState(pDsState);

        ProgramDesc glassDesc;
        glassDesc.addShaderModules(mpScene->getShaderModules());
        glassDesc.addTypeConformances(mpScene->getTypeConformances());
        glassDesc.addShaderLibrary(kShaderFile).vsEntry("vsMain").psEntry("psGlass");
        DefineList glassDefines = mpScene->getSceneDefines();
        glassDefines.add(getShInterpolationDefines());
        mpGlassProgram = Program::create(mpDevice, glassDesc, glassDefines);
        mpGlassVars = ProgramVars::create(mpDevice, mpGlassProgram->getReflector());
        ProgramDesc dumpDesc;
        dumpDesc.addShaderModules(mpScene->getShaderModules());
        dumpDesc.addShaderLibrary(kShaderFile).vsEntry("vsMain").psEntry("psDump");
        mpDumpProgram = Program::create(mpDevice, dumpDesc, staticDefines);
        mpDumpVars = ProgramVars::create(mpDevice, mpDumpProgram->getReflector());
        mpDumpState = GraphicsState::create(mpDevice);
        mpDumpState->setProgram(mpDumpProgram);
        mpDumpState->setRasterizerState(mpRasterState);
        mpDumpState->setDepthStencilState(pDsState);
        // Scene tag of the view dumps: the precompute scene's name, i.e. the loaded scene's file name up to "_withdynamic"
        // (SceneForData_withdynamic_light -> SceneForData, BistroInterior_Wine_withdynamic_lamp -> BistroInterior_Wine).
        // Set on every scene load; edit it in the UI if the precompute scene is named differently.
        {
            const auto paths = mpScene->getImportPaths();
            if (!paths.empty())
            {
                std::string tag = paths.front().stem().string();
                const size_t cut = tag.find("_withdynamic");
                if (cut != std::string::npos && cut > 0) tag = tag.substr(0, cut);
                mViewDumpScene = tag;
            }
        }
        mpGlassState = GraphicsState::create(mpDevice);
        mpGlassState->setProgram(mpGlassProgram);
        mpGlassState->setFbo(mpFbo);
        mpGlassState->setRasterizerState(mpRasterState);
        DepthStencilState::Desc glassDepth;
        glassDepth.setDepthEnabled(true).setDepthWriteMask(true);
        mpGlassState->setDepthStencilState(DepthStencilState::create(glassDepth));
        mpOpaqueResolvePass = ComputePass::create(mpDevice,
            "RenderPasses/AdaptiveSHDemo/ResolveGlassBackground.cs.slang", "main");

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

        mpFilterPass = ComputePass::create(mpDevice, kDynamicFilterFile, "main");

        ProgramDesc compositeDesc;
        compositeDesc.addShaderLibrary(kCompositePassFile)
            .vsEntry("vsMain")
            .psEntry("psMain");

        mpCompositeProgram = Program::create(mpDevice, compositeDesc);
        mpCompositeVars = ProgramVars::create(mpDevice, mpCompositeProgram->getReflector());

        BlendState::Desc blendDesc;
        blendDesc.setRtBlend(0, true)
            .setRtParams(0,
                BlendState::BlendOp::Add,
                BlendState::BlendOp::Add,
                BlendState::BlendFunc::SrcAlpha,
                BlendState::BlendFunc::OneMinusSrcAlpha,
                BlendState::BlendFunc::One,
                BlendState::BlendFunc::OneMinusSrcAlpha);

        DepthStencilState::Desc dsComposite;
        dsComposite.setDepthEnabled(false);

        mpCompositeState = GraphicsState::create(mpDevice);
        mpCompositeState->setProgram(mpCompositeProgram);
        mpCompositeState->setBlendState(BlendState::create(blendDesc));
        mpCompositeState->setDepthStencilState(DepthStencilState::create(dsComposite));
        mpCompositeState->setRasterizerState(mpRasterState);

        mpEmptyVao = Vao::create(Vao::Topology::TriangleStrip);
        mpCompositeState->setVao(mpEmptyVao);

        loadBistroXAtlasLightmapTest(
            mpDevice,
            mpScene
        );
        //REMARK :  set all materials to diffuse for SH testing
        //auto allMat = pScene->getMaterials();

        //for (auto& pMat : allMat)
        //{
        //    // STEP 1: Handle the Base properties (Legacy & Common)
        //    // Since StandardMaterial inherits BasicMaterial, this runs for EVERYONE.
        //    auto pBasicMat = pMat->toBasicMaterial();

        //    if (pBasicMat)
        //    {
        //        // 1. Kill the Specular Color / Shininess
        //        // For Legacy OBJ: This makes it matte.
        //        // For PBR: This ensures the "F0" (Reflectivity at 0 degrees) is black.
        //        pBasicMat->setSpecularParams(float4(0.0f));

        //        // 2. Kill Transmission (Glass/Ghosting)
        //        pBasicMat->setTransmissionColor(float3(0.0f));
        //        pBasicMat->setSpecularTransmission(0.0f);
        //        pBasicMat->setDiffuseTransmission(0.0f);
        //    }

        //    // STEP 2: Handle the PBR-specific properties
        //    // This ONLY runs if the material is actually the modern StandardMaterial type.
        //    StandardMaterial* pStdMat = dynamic_cast<StandardMaterial*>(pMat.get());

        //    if (pStdMat)
        //    {
        //        // 3. Force PBR Roughness (The most important setting for modern renderers)
        //        pStdMat->setRoughness(1.0f);   // 1.0 = Chalk
        //        pStdMat->setMetallic(0.0f);    // 0.0 = Dielectric
        //        pStdMat->setSpecularTransmission(0.0f);
        //        pStdMat->setTransmissionColor(float3(0.0f));
        //    }
        //}
    }
}
