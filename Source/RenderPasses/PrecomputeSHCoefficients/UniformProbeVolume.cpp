#include "UniformProbeVolume.h"
#include <algorithm>
#include <fstream>
#include <iomanip>

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

    float boundsScale = 0.98f;

    // Scale the scene bounds about its center.
    float3 center = 0.5f * (bounds.minPoint + bounds.maxPoint);
    float3 halfExtent = 0.5f * (bounds.maxPoint - bounds.minPoint);
    halfExtent *= boundsScale;

    float3 scaledMin = center - halfExtent;
    float3 scaledMax = center + halfExtent;

    //mMinPoint = bounds.minPoint;
    //mMaxPoint = bounds.maxPoint;
    mMinPoint = scaledMin;
    mMaxPoint = scaledMax;

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

void UniformProbeVolume::saveToFile(const std::string& filename) const
{
    std::ofstream out(filename);
    if (!out)
    {
        logError("Failed to open file for saving: " + filename);
        return;
    }

    out << std::fixed << std::setprecision(8);

    // 1. Header & Grid Configuration
    out << "UNIFORM_GRID_V1\n";

    // Grid Dimensions
    out << mCellResolution.x << " " << mCellResolution.y << " " << mCellResolution.z << "\n";

    // Bounds
    out << mMinPoint.x << " " << mMinPoint.y << " " << mMinPoint.z << "\n";
    out << mMaxPoint.x << " " << mMaxPoint.y << " " << mMaxPoint.z << "\n";
    out << mBuildTimeMs << "\n";
    // 2. Probe Data
    out << "NUM_PROBES " << mProbeData.size() << "\n";

    for (const auto& p : mProbeData)
    {
        // Write Coefficients (9 basis functions * 3 floats)
        // Note: UniformGridCorner stores float4, but we only strictly need float3 for SH.
        // We write .xyz to keep the file smaller/cleaner.
        for (int i = 0; i < 9; ++i)
        {
            out << p.coeffs[i].x << " " << p.coeffs[i].y << " " << p.coeffs[i].z << " ";
        }
        out << "\n";

        // Write Gradients (9 basis functions * 3 channels * 3 dimensions)
        // UniformGridCorner separates gradients by channel (GradR, GradG, GradB)
        for (int i = 0; i < 9; ++i)
        {
            // Grad R (x, y, z)
            out << p.gradR[i].x << " " << p.gradR[i].y << " " << p.gradR[i].z << " ";
            // Grad G (x, y, z)
            out << p.gradG[i].x << " " << p.gradG[i].y << " " << p.gradG[i].z << " ";
            // Grad B (x, y, z)
            out << p.gradB[i].x << " " << p.gradB[i].y << " " << p.gradB[i].z << " ";
        }
        out << "\n";
    }

    out.close();
    logInfo("Successfully saved UniformProbeVolume to " + filename);
}

void UniformProbeVolume::loadFromFile(const std::string& filename)
{
    std::ifstream in(filename);
    if (!in)
    {
        logError("Failed to open file for loading: " + filename);
        return;
    }

    // 1. Check Header
    std::string header;
    in >> header;
    if (header != "UNIFORM_GRID_V1")
    {
        logError("Invalid file format or version mismatch: " + filename);
        return;
    }

    // 2. Load Configuration
    in >> mCellResolution.x >> mCellResolution.y >> mCellResolution.z;
    in >> mMinPoint.x >> mMinPoint.y >> mMinPoint.z;
    in >> mMaxPoint.x >> mMaxPoint.y >> mMaxPoint.z;
    in >> mBuildTimeMs;

    // 3. Re-initialize Internal State (Logic copied from initGrid, but using loaded bounds)
    mProbeCountDim = mCellResolution + uint3(1, 1, 1);
    mTotalProbes = mProbeCountDim.x * mProbeCountDim.y * mProbeCountDim.z;

    float3 dims = float3(std::max(1u, mCellResolution.x), std::max(1u, mCellResolution.y), std::max(1u, mCellResolution.z));
    mCellSize = (mMaxPoint - mMinPoint) / dims;

    // Resize Storage
    mProbeData.clear();
    mProbeData.resize(mTotalProbes);

    // 4. Load Probe Data
    std::string tag;
    size_t numProbesInFile;
    in >> tag >> numProbesInFile;

    if (numProbesInFile != mTotalProbes)
    {
        logWarning("Probe count in file (" + std::to_string(numProbesInFile) +
            ") does not match calculated grid size (" + std::to_string(mTotalProbes) + "). File might be corrupt or configuration mismatched.");
    }

    // Ensure we don't overflow if file has more, or underflow if less
    size_t probesToRead = std::min((size_t)mTotalProbes, numProbesInFile);

    for (size_t k = 0; k < probesToRead; ++k)
    {
        UniformGridCorner& p = mProbeData[k];

        // Read Coefficients
        for (int i = 0; i < 9; ++i)
        {
            float3 val;
            in >> val.x >> val.y >> val.z;
            p.coeffs[i] = float4(val, 0.0f); // Restore w=0 padding
        }

        // Read Gradients
        for (int i = 0; i < 9; ++i)
        {
            float3 gR, gG, gB;
            in >> gR.x >> gR.y >> gR.z;
            in >> gG.x >> gG.y >> gG.z;
            in >> gB.x >> gB.y >> gB.z;

            p.gradR[i] = float4(gR, 0.0f);
            p.gradG[i] = float4(gG, 0.0f);
            p.gradB[i] = float4(gB, 0.0f);
        }
    }

    in.close();
    logInfo("Successfully loaded UniformProbeVolume from " + filename);
}
