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
#include <fstream>
#include "AdaptiveSHDemo.h"
#include "Rendering/Lights/EmissivePowerSampler.h"
#include "Rendering/Lights/EmissiveUniformSampler.h"
#include "Utils/Logger.h"

#define PROBE_MODE_ADAPTIVE 0
#define PROBE_MODE_UNIFORM  1
 // CHANGE THIS LINE TO SWITCH MODES:
//#define CURRENT_PROBE_MODE PROBE_MODE_UNIFORM
#define CURRENT_PROBE_MODE PROBE_MODE_ADAPTIVE

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
const std::string loadFromFileName = "DirectAbsErr2N6EdgeGradientMetricDataScene4096spp.txt";
//const std::string loadFromFileName = "DirectAbsErr2EdgeMetricCornellThinSlabV2.txt";
//const std::string loadFromFileName = "DirectAbsErr2HessianMetricCornellThinSlabV2.txt";
//const std::string loadFromFileName = "DirectAbsErr2N6HessianMetricDataScene4096spp.txt";
const char kShaderFile[] = "RenderPasses/AdaptiveSHDemo/AdaptiveGridShader.slang";

#else
//const std::string loadFromFileName = "U32DataScene.txt";
//const std::string loadFromFileName = "U32DataScene_4096spp_4spt.txt";
//const std::string loadFromFileName = "U64CornellShadowBoundaryScene.txt";
//const std::string loadFromFileName = "U32CornellShadowBoundaryScene.txt";
//const std::string loadFromFileName = "U64DataScene.txt";
const std::string loadFromFileName = "U64DataScene_4096spp.txt";
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
    mpFbo = Fbo::create(mpDevice);
    Sampler::Desc samplerDesc;
    samplerDesc.setFilterMode(TextureFilteringMode::Linear, TextureFilteringMode::Linear, TextureFilteringMode::Linear);
    mpLinearSampler = mpDevice->createSampler(samplerDesc);
}

namespace
{
    constexpr float kPi = 3.14159265358979323846f;

    float degreesToRadians(float degrees)
    {
        return degrees * kPi / 180.0f;
    }

    float radiansToDegrees(float radians)
    {
        return radians * 180.0f / kPi;
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
    static uint32_t frameIndex = 0;
    static constexpr uint32_t kWarmupFrames = 100;
    static constexpr uint32_t kMeasureFrames = 1000;

    static auto lastTime = std::chrono::high_resolution_clock::now();
    static std::vector<double> frameTimesMs;

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
        // ------------------------------------------------------------------
        // PASS 1: STATIC GEOMETRY (lightmaps)
        // ------------------------------------------------------------------
        auto applyVar = mpStaticVars->getRootVar();
        bindDataSceneData(applyVar);
        //bindCornellData(applyVar);
        //bindCornellVisibilitySlabData(applyVar);
#if CURRENT_PROBE_MODE == PROBE_MODE_ADAPTIVE
        applyVar["gCornerBuffer"] = mAdaptiveProbeVolume->getCornerBuffer();
        applyVar["gProbeBuffer"] = mAdaptiveProbeVolume->getProbeBuffer();

        applyVar["gSeedProbeIndices"] = mAdaptiveProbeVolume->getSeedProbeIndexBuffer();

        applyVar["SeedGridCB"]["gUseSeedGrid"] = mAdaptiveProbeVolume->getUseSeedGrid() ? 1u : 0u;
        applyVar["SeedGridCB"]["gSeedMinPoint"] = mAdaptiveProbeVolume->getSeedMinPoint();
        applyVar["SeedGridCB"]["gSeedCellSize"] = mAdaptiveProbeVolume->getSeedCellSize();
        applyVar["SeedGridCB"]["gSeedResolution"] = mAdaptiveProbeVolume->getSeedResolution();

#else
        mUniformProbeVolume->bindShaderData(applyVar);
#endif

        mpScene->rasterize(pRenderContext, mpGraphicsState.get(), mpStaticVars.get(), mpRasterState, mpRasterState);

        auto now = std::chrono::high_resolution_clock::now();
        double frameMs = std::chrono::duration<double, std::milli>(now - lastTime).count();
        lastTime = now;

        if (frameIndex >= kWarmupFrames && frameIndex < kWarmupFrames + kMeasureFrames)
        {
            frameTimesMs.push_back(frameMs);
        }

        if (frameIndex == kWarmupFrames + kMeasureFrames)
        {
            double sum = 0.0;
            for (double v : frameTimesMs) sum += v;

            double meanMs = sum / double(frameTimesMs.size());
            double meanFps = 1000.0 / meanMs;

            double variance = 0.0;
            for (double v : frameTimesMs)
            {
                double d = v - meanMs;
                variance += d * d;
            }
            variance /= double(frameTimesMs.size());
            double stdMs = std::sqrt(variance);
            uint64_t probeCount = 0;
#if CURRENT_PROBE_MODE == PROBE_MODE_ADAPTIVE
            std::string modeName = "Adaptive";
            probeCount = mAdaptiveProbeVolume->getProbes().size();
#else
            std::string modeName = "Uniform";
            const uint3 res = mUniformProbeVolume->getCellResolution();
            probeCount = uint64_t(res.x) * uint64_t(res.y) * uint64_t(res.z);
#endif

            //std::ofstream out("AdaptiveSHDemo_RuntimeFPS.csv", std::ios::app);
            //out << modeName << ","
            //    << loadFromFileName << ","
            //    << probeCount << ","
            //    << kWarmupFrames << ","
            //    << kMeasureFrames << ","
            //    << std::fixed << std::setprecision(4)
            //    << meanMs << ","
            //    << stdMs << ","
            //    << meanFps << "\n";

            //out.close();

            //logInfo("Runtime measurement finished: file={}, mean {:.4f} ms, std {:.4f} ms, FPS {:.2f}",
            //    loadFromFileName, meanMs, stdMs, meanFps);
        }

        frameIndex++;

        if (mbShowAdaptiveGrid)
        {
            mpProbeVisualizePass->setCameraData(
                mpScene->getCamera()->getViewProjMatrix()
            );

            mpProbeVisualizePass->execute(pRenderContext, mpFbo);
        }
    }
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

void AdaptiveSHDemo::renderUI(Gui::Widgets& widget) {
    widget.text("Loaded probe file: " + loadFromFileName);

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

void AdaptiveSHDemo::loadLightmaps()
{
    // Load as a 2D texture. Falcor handles EXR (HDR) automatically.
    // We set loadAsSrgb to false because lightmaps contain linear radiance data.
    //auto loadOne = [&](size_t idx, ref<Texture>& dst, const std::string& debugName)
    //    {
    //        if (idx >= mBakeTargets.size())
    //        {
    //            logWarning("Bake target index {} is out of range.", idx);
    //            return;
    //        }

    //        const auto& target = mBakeTargets[idx];

    //        dst = Texture::createFromFile(
    //            mpDevice,
    //            target.outputPath,
    //            true,
    //            false,
    //            ResourceBindFlags::ShaderResource
    //        );

    //        if (dst)
    //        {
    //            dst->setName(debugName);
    //            logInfo("Successfully loaded lightmap '{}' from {}", target.name, target.outputPath);
    //        }
    //        else
    //        {
    //            logWarning("Failed to load lightmap '{}' from {}", target.name, target.outputPath);
    //        }
    //    };

    //loadOne(0, mpFloorLightmap, "FloorLightmap");
    //loadOne(1, mpLeftWallLightmap, "LeftWallLightmap");
    //loadOne(2, mpRightWallLightmap, "RightWallLightmap");
    //loadOne(3, mpRoofLeftLightmap, "RoofLeftLightmap");
    //loadOne(4, mpRoofRightLightmap, "RoofRightLightmap");

    //loadOne(5, mpPillar0Lightmap, "Pillar0Lightmap");
    //loadOne(6, mpPillar1Lightmap, "Pillar1Lightmap");
    //loadOne(7, mpPillar2Lightmap, "Pillar2Lightmap");
    //loadOne(8, mpPillar3Lightmap, "Pillar3Lightmap");
    //loadOne(9, mpPillar4Lightmap, "Pillar4Lightmap");
    //loadOne(10, mpPillar5Lightmap, "Pillar5Lightmap");
    //loadOne(11, mpPillar6Lightmap, "Pillar6Lightmap");
    //loadOne(12, mpPillar7Lightmap, "Pillar7Lightmap");
}

void AdaptiveSHDemo::setScene(RenderContext* pRenderContext, const ref<Scene>& pScene)
{
    mpScene = pScene;
    if (mpScene)
    {
        std::ifstream check("AdaptiveSHDemo_RuntimeFPS.csv");
        if (!check.good())
        {
            std::ofstream out("AdaptiveSHDemo_RuntimeFPS.csv");
            out << "Mode,ProbeFile,ProbeCount,WarmupFrames,MeasuredFrames,MeanFrameMs,StdFrameMs,MeanFPS\n";
        }

        setupDataSceneBakeTargets();
        //setupCornellBakeTargets();
        //setupCornellVisibilitySlabBakeTargets();
        //init probe visual pass
        mpProbeVisualizePass = ProbeVisualizePass::create(mpDevice, mpScene->getSceneDefines());
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
        mAdaptiveProbeVolume = AdaptiveProbeVolume::create(mpDevice);

            mAdaptiveProbeVolume->loadFromFile(loadFromFileName);
            mAdaptiveProbeVolume->uploadToGPU();
            mpProbeVisualizePass->setVolumeData(mAdaptiveProbeVolume->getProbes());
#else
        mUniformProbeVolume = UniformProbeVolume::create(mpDevice);
            mUniformProbeVolume->loadFromFile(loadFromFileName);
            mUniformProbeVolume->uploadToGPU();
            mpProbeVisualizePass->setUniformVolumeData(
                mUniformProbeVolume->getMinPoint(),
                mUniformProbeVolume->getCellSize(),
                mUniformProbeVolume->getCellResolution()
            );
#endif
      

        ProgramDesc previewDesc;
        previewDesc.addShaderModules(mpScene->getShaderModules());
        previewDesc.addShaderLibrary(kShaderFile)
            .vsEntry("vsMain")
            .psEntry("psMain");

        mpStaticProgram = Program::create(mpDevice, previewDesc, mpScene->getSceneDefines());
        mpStaticVars = ProgramVars::create(mpDevice, mpStaticProgram->getReflector());

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

        //loadLightmaps();
        loadDataSceneLightmaps();
        //loadCornellLightmaps();
        //loadCornellVisibilitySlabLightmaps();
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

void AdaptiveSHDemo::loadDataSceneLightmaps()
{
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
                logInfo("Successfully loaded data-scene lightmap '{}' from {}", target.name, target.outputPath);
            }
            else
            {
                logWarning("Failed to load data-scene lightmap '{}' from {}", target.name, target.outputPath);
            }
        };

    loadOne(0, mpDataFloorLightmap, "DataFloorLightmap");
    loadOne(1, mpDataCeilingLightmap, "DataCeilingLightmap");
    loadOne(2, mpDataLeftWallLightmap, "DataLeftWallLightmap");
    loadOne(3, mpDataRightWallLightmap, "DataRightWallLightmap");
    loadOne(4, mpDataBackWallLightmap, "DataBackWallLightmap");
    loadOne(5, mpDataFrontWallLightmap, "DataFrontWallLightmap");

    loadOne(6, mpDataTallBoxALightmap, "DataTallBoxALightmap");
    loadOne(7, mpDataWideBoxBLightmap, "DataWideBoxBLightmap");
    loadOne(8, mpDataBlockCLightmap, "DataBlockCLightmap");
    loadOne(9, mpDataLowBoxDLightmap, "DataLowBoxDLightmap");
    loadOne(10, mpDataThinSlabELightmap, "DataThinSlabELightmap");
    loadOne(11, mpDataThinSlabFLightmap, "DataThinSlabFLightmap");
    loadOne(12, mpDataTallBoxJLightmap, "DataTallBoxJLightmap");
    loadOne(13, mpDataShortBoxLLightmap, "DataShortBoxLLightmap");
}

void AdaptiveSHDemo::bindDataSceneData(ShaderVar applyVar)
{
    applyVar["gLinearSampler"] = mpLinearSampler;

    applyVar["gFloorLightmap"] = mpDataFloorLightmap;
    applyVar["gCeilingLightmap"] = mpDataCeilingLightmap;
    applyVar["gLeftWallLightmap"] = mpDataLeftWallLightmap;
    applyVar["gRightWallLightmap"] = mpDataRightWallLightmap;
    applyVar["gBackWallLightmap"] = mpDataBackWallLightmap;
    applyVar["gFrontWallLightmap"] = mpDataFrontWallLightmap;

    applyVar["gTallBoxALightmap"] = mpDataTallBoxALightmap;
    applyVar["gWideBoxBLightmap"] = mpDataWideBoxBLightmap;
    applyVar["gBlockCLightmap"] = mpDataBlockCLightmap;
    applyVar["gLowBoxDLightmap"] = mpDataLowBoxDLightmap;
    applyVar["gThinSlabELightmap"] = mpDataThinSlabELightmap;
    applyVar["gThinSlabFLightmap"] = mpDataThinSlabFLightmap;
    applyVar["gTallBoxJLightmap"] = mpDataTallBoxJLightmap;
    applyVar["gShortBoxLLightmap"] = mpDataShortBoxLLightmap;

    applyVar["PerFrameCB"]["gFloorInstanceID"] = 0;
    applyVar["PerFrameCB"]["gCeilingInstanceID"] = 1;
    applyVar["PerFrameCB"]["gLeftWallInstanceID"] = 2;
    applyVar["PerFrameCB"]["gRightWallInstanceID"] = 3;
    applyVar["PerFrameCB"]["gBackWallInstanceID"] = 4;
    applyVar["PerFrameCB"]["gFrontWallInstanceID"] = 5;

    applyVar["PerFrameCB"]["gTallBoxAInstanceID"] = 6;
    applyVar["PerFrameCB"]["gWideBoxBInstanceID"] = 7;
    applyVar["PerFrameCB"]["gBlockCInstanceID"] = 8;
    applyVar["PerFrameCB"]["gLowBoxDInstanceID"] = 9;
    applyVar["PerFrameCB"]["gThinSlabEInstanceID"] = 10;
    applyVar["PerFrameCB"]["gThinSlabFInstanceID"] = 11;
    applyVar["PerFrameCB"]["gTallBoxJInstanceID"] = 12;
    applyVar["PerFrameCB"]["gShortBoxLInstanceID"] = 13;

    applyVar["PerFrameCB"]["gTallBoxACenterW"] = float3(-3.8f, 1.3f, -2.8f);
    applyVar["PerFrameCB"]["gWideBoxBCenterW"] = float3(-1.9f, 0.5f, 2.2f);
    applyVar["PerFrameCB"]["gBlockCCenterW"] = float3(2.8f, 0.8f, -2.6f);
    applyVar["PerFrameCB"]["gLowBoxDCenterW"] = float3(1.4f, 0.35f, 2.8f);
    applyVar["PerFrameCB"]["gThinSlabECenterW"] = float3(-0.4f, 1.5f, -0.8f);
    applyVar["PerFrameCB"]["gThinSlabFCenterW"] = float3(3.4f, 1.2f, 0.8f);
    applyVar["PerFrameCB"]["gTallBoxJCenterW"] = float3(-4.2f, 1.5f, 2.8f);
    applyVar["PerFrameCB"]["gShortBoxLCenterW"] = float3(-0.8f, 0.45f, -3.3f);

    applyVar["PerFrameCB"]["gTallBoxAHalfExtentW"] = float3(0.45f, 1.30f, 0.45f);
    applyVar["PerFrameCB"]["gWideBoxBHalfExtentW"] = float3(0.90f, 0.50f, 0.60f);
    applyVar["PerFrameCB"]["gBlockCHalfExtentW"] = float3(0.72f, 0.80f, 0.72f);
    applyVar["PerFrameCB"]["gLowBoxDHalfExtentW"] = float3(1.05f, 0.35f, 1.00f);
    applyVar["PerFrameCB"]["gThinSlabEHalfExtentW"] = float3(0.175f, 1.50f, 0.70f);
    applyVar["PerFrameCB"]["gThinSlabFHalfExtentW"] = float3(0.36f, 1.20f, 0.85f);
    applyVar["PerFrameCB"]["gTallBoxJHalfExtentW"] = float3(0.75f, 1.25f, 0.75f);
    applyVar["PerFrameCB"]["gShortBoxLHalfExtentW"] = float3(0.75f, 0.45f, 0.75f);
}

void AdaptiveSHDemo::setupDataSceneBakeTargets()
{
    mBakeTargets =
    {
        // Room shell.
        //{ "Floor",     0, 1024, 1024, "BakedData_Floor.exr"     },
        { "Floor",     0, 2048, 2048, "BakedData_Floor.exr"     },
        { "Ceiling",   1, 1024, 1024, "BakedData_Ceiling.exr"   },
        { "LeftWall",  2, 1024,  512, "BakedData_LeftWall.exr"  },
        { "RightWall", 3, 1024,  512, "BakedData_RightWall.exr" },
        { "BackWall",  4, 1024,  512, "BakedData_BackWall.exr"  },
        { "FrontWall", 5, 1024,  512, "BakedData_FrontWall.exr" },

        //Box targets.
        //halfExtent = scaling * 0.5
        //rotationEulerDeg = same as SceneForData.pyscene

       { "TallBoxA", 6, 512, 512, "BakedData_TallBoxA.exr",
           BakeTargetType::Pillar,
           float3(-3.8f, 1.3f, -2.8f),
           float3(0.45f, 1.30f, 0.45f),
           float3(0.0f, 0.0f, 0.0f) },

       { "WideBoxB", 7, 1024, 1024, "BakedData_WideBoxB.exr",
           BakeTargetType::Pillar,
           float3(-1.9f, 0.5f, 2.2f),
           float3(0.90f, 0.50f, 0.60f),
           float3(0.0f, 0.0f, 0.0f) },

       { "BlockC", 8, 512, 512, "BakedData_BlockC.exr",
           BakeTargetType::Pillar,
           float3(2.8f, 0.8f, -2.6f),
           float3(0.55f, 0.80f, 0.55f),
           float3(0.0f, 18.0f, 0.0f) },

       { "LowBoxD", 9, 512, 512, "BakedData_LowBoxD.exr",
           BakeTargetType::Pillar,
           float3(1.4f, 0.35f, 2.8f),
           float3(0.80f, 0.35f, 0.70f),
           float3(0.0f, -22.0f, 0.0f) },

       { "ThinSlabE", 10, 512, 512, "BakedData_ThinSlabE.exr",
           BakeTargetType::Pillar,
           float3(-0.4f, 1.5f, -0.8f),
           float3(0.175f, 1.50f, 0.70f),
           float3(0.0f, 0.0f, 0.0f) },

       { "ThinSlabF", 11, 512, 512, "BakedData_ThinSlabF.exr",
           BakeTargetType::Pillar,
           float3(3.4f, 1.2f, 0.8f),
           float3(0.175f, 1.20f, 0.80f),
           float3(0.0f, 12.0f, 0.0f) },

       { "TallBoxJ", 12, 512, 512, "BakedData_TallBoxJ.exr",
           BakeTargetType::Pillar,
           float3(-4.2f, 1.5f, 2.8f),
           float3(0.35f, 1.10f, 0.35f),
           float3(0.0f, -18.0f, 8.0f) },

       { "ShortBoxL", 13, 512, 512, "BakedData_ShortBoxL.exr",
           BakeTargetType::Pillar,
           float3(-0.8f, 0.45f, -3.3f),
           float3(0.50f, 0.45f, 0.50f),
           float3(0.0f, 40.0f, 0.0f) },
    };
}

void AdaptiveSHDemo::setupCornellBakeTargets()
{
    mBakeTargets =
    {
        // Room shell targets (Quads)
        { "Floor",     0, 1024, 1024, "BakedCornell_Floor.exr"     },
        { "Ceiling",   1, 1024, 1024, "BakedCornell_Ceiling.exr"   },
        { "BackWall",  2, 1024, 1024, "BakedCornell_BackWall.exr"  },
        { "LeftWall",  3, 1024, 1024, "BakedCornell_LeftWall.exr"  },
        { "RightWall", 4, 1024, 1024, "BakedCornell_RightWall.exr" },

        // Light is instance ID 5; we skip baking the emissive light source.

        // Visibility-discontinuity test slab (Cube/Pillar)
        { "VisibilitySlab", 6, 512, 512, "BakedCornell_VisibilitySlab.exr",
            BakeTargetType::Pillar,
            float3(0.12f, 0.28f, 0.00f),       // translation
            float3(0.15f, 0.0175f, 0.11f),     // half extent (scaling * 0.5)
            float3(0.0f, 0.0f, 0.0f)           // rotation
        }
    };
}

void AdaptiveSHDemo::loadCornellLightmaps()
{
    auto loadOne = [&](size_t idx, ref<Texture>& dst, const std::string& debugName)
        {
            if (idx >= mBakeTargets.size()) return;
            const auto& target = mBakeTargets[idx];
            dst = Texture::createFromFile(mpDevice, target.outputPath, true, false, ResourceBindFlags::ShaderResource);
            if (dst) dst->setName(debugName);
        };

    loadOne(0, mpCornellFloorLightmapShadowBoundaryTestScene, "CornellFloorLightmap");
    loadOne(1, mpCornellCeilingLightmapShadowBoundaryTestScene, "CornellCeilingLightmap");
    loadOne(2, mpCornellBackWallLightmapShadowBoundaryTestScene, "CornellBackWallLightmap");
    loadOne(3, mpCornellLeftWallLightmapShadowBoundaryTestScene, "CornellLeftWallLightmap");
    loadOne(4, mpCornellRightWallLightmapShadowBoundaryTestScene, "CornellRightWallLightmap");
    loadOne(5, mpCornellThinSlabLightmapShadowBoundaryTestScene, "CornellVisibilitySlabLightmap");
}

void AdaptiveSHDemo::bindCornellData(ShaderVar applyVar)
{
    applyVar["gLinearSampler"] = mpLinearSampler;

    // Bind Textures
    applyVar["gFloorLightmap"] = mpCornellFloorLightmapShadowBoundaryTestScene;
    applyVar["gCeilingLightmap"] = mpCornellCeilingLightmapShadowBoundaryTestScene;
    applyVar["gBackWallLightmap"] = mpCornellBackWallLightmapShadowBoundaryTestScene;
    applyVar["gLeftWallLightmap"] = mpCornellLeftWallLightmapShadowBoundaryTestScene;
    applyVar["gRightWallLightmap"] = mpCornellRightWallLightmapShadowBoundaryTestScene;
    applyVar["gVisibilitySlabLightmap"] = mpCornellThinSlabLightmapShadowBoundaryTestScene;

    // Bind Instance IDs
    applyVar["PerFrameCB"]["gFloorInstanceID"] = 0;
    applyVar["PerFrameCB"]["gCeilingInstanceID"] = 1;
    applyVar["PerFrameCB"]["gBackWallInstanceID"] = 2;
    applyVar["PerFrameCB"]["gLeftWallInstanceID"] = 3;
    applyVar["PerFrameCB"]["gRightWallInstanceID"] = 4;
    applyVar["PerFrameCB"]["gVisibilitySlabInstanceID"] = 6;

    // Bind Pillar Data for Slab
    applyVar["PerFrameCB"]["gVisibilitySlabCenterW"] = float3(0.12f, 0.28f, 0.00f);
    applyVar["PerFrameCB"]["gVisibilitySlabHalfExtentW"] = float3(0.15f, 0.0175f, 0.11f);
}

void AdaptiveSHDemo::setupCornellVisibilitySlabBakeTargets()
{
    mBakeTargets =
    {
        // Room shell targets (Quads)
        { "Floor",     0, 1024, 1024, "BakedCornellSlab_Floor.exr"     },
        { "Ceiling",   1, 1024, 1024, "BakedCornellSlab_Ceiling.exr"   },
        { "BackWall",  2, 1024, 1024, "BakedCornellSlab_BackWall.exr"  },
        { "LeftWall",  3, 1024, 1024, "BakedCornellSlab_LeftWall.exr"  },
        { "RightWall", 4, 1024, 1024, "BakedCornellSlab_RightWall.exr" },

        // Visibility-discontinuity test slab (Cube/Pillar)
        { "VisibilitySlab", 6, 512, 512, "BakedCornellSlab_VisibilitySlab.exr",
            BakeTargetType::Pillar,
            float3(0.0f, 0.21f, -0.02f),       // Center (slabX, slabY, slabZ)
            float3(0.01f, 0.21f, 0.18f),       // Half extent (scaling * 0.5)
            float3(0.0f, 0.0f, 0.0f)           // Rotation
        }
    };
}

void AdaptiveSHDemo::loadCornellVisibilitySlabLightmaps()
{
    auto loadOne = [&](size_t idx, ref<Texture>& dst, const std::string& debugName)
        {
            if (idx >= mBakeTargets.size()) return;
            const auto& target = mBakeTargets[idx];
            dst = Texture::createFromFile(mpDevice, target.outputPath, true, false, ResourceBindFlags::ShaderResource);
            if (dst) dst->setName(debugName);
        };

    loadOne(0, mpCornellFloorLightmapVisibilitySlab, "CornellFloorLightmapSlab");
    loadOne(1, mpCornellCeilingLightmapVisibilitySlab, "CornellCeilingLightmapSlab");
    loadOne(2, mpCornellBackWallLightmapVisibilitySlab, "CornellBackWallLightmapSlab");
    loadOne(3, mpCornellLeftWallLightmapVisibilitySlab, "CornellLeftWallLightmapSlab");
    loadOne(4, mpCornellRightWallLightmapVisibilitySlab, "CornellRightWallLightmapSlab");
    loadOne(5, mpCornellSlabLightmapVisibilitySlab, "CornellVisibilitySlabLightmapSlab");
}

void AdaptiveSHDemo::bindCornellVisibilitySlabData(ShaderVar applyVar)
{
    applyVar["gLinearSampler"] = mpLinearSampler;

    // Bind Textures
    applyVar["gFloorLightmap"] = mpCornellFloorLightmapVisibilitySlab;
    applyVar["gCeilingLightmap"] = mpCornellCeilingLightmapVisibilitySlab;
    applyVar["gBackWallLightmap"] = mpCornellBackWallLightmapVisibilitySlab;
    applyVar["gLeftWallLightmap"] = mpCornellLeftWallLightmapVisibilitySlab;
    applyVar["gRightWallLightmap"] = mpCornellRightWallLightmapVisibilitySlab;
    applyVar["gVisibilitySlabLightmap"] = mpCornellSlabLightmapVisibilitySlab;

    // Bind Instance IDs
    applyVar["PerFrameCB"]["gFloorInstanceID"] = 0;
    applyVar["PerFrameCB"]["gCeilingInstanceID"] = 1;
    applyVar["PerFrameCB"]["gBackWallInstanceID"] = 2;
    applyVar["PerFrameCB"]["gLeftWallInstanceID"] = 3;
    applyVar["PerFrameCB"]["gRightWallInstanceID"] = 4;
    applyVar["PerFrameCB"]["gVisibilitySlabInstanceID"] = 6;

    // Bind Pillar Data for Slab
    applyVar["PerFrameCB"]["gVisibilitySlabCenterW"] = float3(0.0f, 0.21f, -0.02f);
    applyVar["PerFrameCB"]["gVisibilitySlabHalfExtentW"] = float3(0.01f, 0.21f, 0.18f);
}
