#include "UniformProbeVolume.h"
#include <algorithm>

ref<UniformProbeVolume> UniformProbeVolume::create(ref<Device> pDevice)
{
    return ref<UniformProbeVolume>(new UniformProbeVolume(pDevice));
}

UniformProbeVolume::UniformProbeVolume(ref<Device> pDevice)
    : mpDevice(pDevice)
{
}

void UniformProbeVolume::initGrid(const ref<Scene>& pScene, uint3 cellResolution)
{
    mCellResolution = cellResolution;

    // CORRECTION: Probes are at corners, so we need (N + 1)
    mProbeCountDim = mCellResolution + uint3(1, 1, 1);
    mTotalProbes = mProbeCountDim.x * mProbeCountDim.y * mProbeCountDim.z;

    auto bounds = pScene->getSceneBounds();
    mMinPoint = bounds.minPoint;
    mMaxPoint = bounds.maxPoint;

    // Cell size is based on the volume bounds divided by CELL count
    float3 dims = float3(std::max(1u, mCellResolution.x), std::max(1u, mCellResolution.y), std::max(1u, mCellResolution.z));
    mCellSize = (mMaxPoint - mMinPoint) / dims;

    // Resize Storage
    mProbeData.clear();
    mProbeData.resize(mTotalProbes);
    std::memset(mProbeData.data(), 0, mTotalProbes * sizeof(UniformGridCorner));
}

void UniformProbeVolume::getProbePositions(std::vector<float3>& outPositions) const
{
    outPositions.clear();
    outPositions.reserve(mTotalProbes);

    // Generate positions for every CORNER
    // Z -> Y -> X iteration order
    for (uint z = 0; z < mProbeCountDim.z; ++z)
    {
        for (uint y = 0; y < mProbeCountDim.y; ++y)
        {
            for (uint x = 0; x < mProbeCountDim.x; ++x)
            {
                // Position = Origin + Index * CellSize
                float3 pos = mMinPoint + float3(x, y, z) * mCellSize;
                outPositions.push_back(pos);
            }
        }
    }
}

void UniformProbeVolume::setProbeData(uint32_t probeIndex, const std::vector<float3>& coeffs, const std::vector<GradSHCoeff>& grads)
{
    if (probeIndex >= mTotalProbes) return;

    UniformGridCorner& p = mProbeData[probeIndex];
    size_t count = std::min((size_t)9, coeffs.size());

    for (size_t i = 0; i < 9; ++i)
    {
        if (i < count)
        {
            p.coeffs[i] = float4(coeffs[i], 0.0f);
            p.gradR[i] = float4(grads[i].r, 0.0f);
            p.gradG[i] = float4(grads[i].g, 0.0f);
            p.gradB[i] = float4(grads[i].b, 0.0f);
        }
        else
        {
            p.coeffs[i] = float4(0.f);
            p.gradR[i] = float4(0.f); p.gradG[i] = float4(0.f); p.gradB[i] = float4(0.f);
        }
    }
}

void UniformProbeVolume::uploadToGPU()
{
    if (mTotalProbes == 0) return;

    mpProbeDataBuffer = mpDevice->createStructuredBuffer(
        sizeof(UniformGridCorner),
        (uint32_t)mProbeData.size(),
        ResourceBindFlags::ShaderResource,
        MemoryType::DeviceLocal,
        mProbeData.data()
    );
    mpProbeDataBuffer->setName("UniformProbeBuffer");
}

void UniformProbeVolume::bindShaderData(ShaderVar& var)
{
    var["gUniformProbes"] = mpProbeDataBuffer;
    var["PerFrameCB"]["gGridMin"] = mMinPoint;
    //var["PerFrameCB"]["gGridMax"] = mMaxPoint;
    // Pass the PROBE DIMENSION (Nodes), not the Cell Dimension, for stride calculations
    var["PerFrameCB"]["gProbeCountDim"] = mProbeCountDim;
    var["PerFrameCB"]["gCellSize"] = mCellSize;
}
