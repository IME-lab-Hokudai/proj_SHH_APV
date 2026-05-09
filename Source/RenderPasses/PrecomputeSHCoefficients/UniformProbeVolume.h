#pragma once
#include "Falcor.h"
#include "ProbeSamplingData.slang" // Defines UniformGridCorner
#include "envMap_SH.h"

using namespace Falcor;

class UniformProbeVolume : public Object
{
public:
    static ref<UniformProbeVolume> create(ref<Device> pDevice);

    // 1. Initialization
    // resolution = Number of CELLS (e.g. 16x16x16)
    // The actual number of probes will be (17x17x17)
    void initGrid(const ref<Scene>& pScene, uint3 cellResolution);

    // 2. Baking Interface
    void getProbePositions(std::vector<float3>& outPositions) const;
    void setProbeData(uint32_t probeIndex, const std::vector<float3>& coeffs, const std::vector<GradSHCoeff>& grads);

    // 3. GPU Management
    void uploadToGPU();
    void bindShaderData(ShaderVar& var);

    // Getters
    uint3 getCellResolution() const { return mProbeResolution; } // 16,16,16
    uint3 getProbeCountDim() const { return mCornerResolution; }   // 17,17,17
    float3 getMinPoint() const { return mMinPoint; }
    float3 getMaxPoint() const { return mMaxPoint; }
    float3 getCellSize() const { return mCellSize; }

    void saveToFile(const std::string& filename) const;
    void loadFromFile(const std::string& filename);

    double getBuildTimeMs() const { return mBuildTimeMs; }
    void setBuildTimeMs(double t) { mBuildTimeMs = t; }

private:
    UniformProbeVolume(ref<Device> pDevice);

    ref<Device> mpDevice;

    // Dimensions
    uint3 mProbeResolution = { 0, 0, 0 }; // Number of Boxes
    uint3 mCornerResolution = { 0, 0, 0 };  // Number of Corners (Res + 1)
    uint32_t mTotalProbes = 0;           // Flattened count

    float3 mMinPoint = { 0, 0, 0 };
    float3 mMaxPoint = { 0, 0, 0 };
    float3 mCellSize = { 0, 0, 0 };

    // Unified Storage
    std::vector<UniformGridCorner> mProbeData;
    ref<Buffer> mpProbeDataBuffer;

    //statistic
    double mBuildTimeMs = 0.0;
};
