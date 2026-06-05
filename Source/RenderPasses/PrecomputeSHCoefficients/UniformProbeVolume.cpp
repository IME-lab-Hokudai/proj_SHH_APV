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
    mProbeResolution = cellResolution;

    // CORRECTION: Probes are at corners, so we need (N + 1)
    mCornerResolution = mProbeResolution + uint3(1, 1, 1);
    mTotalProbes = mCornerResolution.x * mCornerResolution.y * mCornerResolution.z;

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
    float3 dims = float3(std::max(1u, mProbeResolution.x), std::max(1u, mProbeResolution.y), std::max(1u, mProbeResolution.z));
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
    for (uint z = 0; z < mCornerResolution.z; ++z)
    {
        for (uint y = 0; y < mCornerResolution.y; ++y)
        {
            for (uint x = 0; x < mCornerResolution.x; ++x)
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
    var["PerFrameCB"]["gProbeCountDim"] = mCornerResolution;
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
    out << mProbeResolution.x << " " << mProbeResolution.y << " " << mProbeResolution.z << "\n";

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
    in >> mProbeResolution.x >> mProbeResolution.y >> mProbeResolution.z;
    in >> mMinPoint.x >> mMinPoint.y >> mMinPoint.z;
    in >> mMaxPoint.x >> mMaxPoint.y >> mMaxPoint.z;
    in >> mBuildTimeMs;
    // 3. Re-initialize Internal State (Logic copied from initGrid, but using loaded bounds)
    mCornerResolution = mProbeResolution + uint3(1, 1, 1);
    mTotalProbes = mCornerResolution.x * mCornerResolution.y * mCornerResolution.z;

    float3 dims = float3(std::max(1u, mProbeResolution.x), std::max(1u, mProbeResolution.y), std::max(1u, mProbeResolution.z));
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


namespace
{
    static float clamp01CPU(float v)
    {
        return std::max(0.0f, std::min(1.0f, v));
    }

    static float hermite1D_CPU(float t, float v0, float s0, float v1, float s1)
    {
        float t2 = t * t;
        float t3 = t2 * t;

        float h00 = 2.0f * t3 - 3.0f * t2 + 1.0f;
        float h10 = t3 - 2.0f * t2 + t;
        float h01 = -2.0f * t3 + 3.0f * t2;
        float h11 = t3 - t2;

        return h00 * v0 + h10 * s0 + h01 * v1 + h11 * s1;
    }

    static float3 hermite1D_CPU(float t, float3 v0, float3 s0, float3 v1, float3 s1)
    {
        return float3(
            hermite1D_CPU(t, v0.x, s0.x, v1.x, s1.x),
            hermite1D_CPU(t, v0.y, s0.y, v1.y, s1.y),
            hermite1D_CPU(t, v0.z, s0.z, v1.z, s1.z)
        );
    }

    static float3 evaluateIrradiance(float3 n, const std::vector<float3>& shCoeffs)
    {
        // Real Spherical Harmonics up to band 2
        // Coordinate convention: x = forward, y = right, z = up
        // so we need to remap falcor coords (Y-up) to polar coords (Z-up)
        float x = n.z;
        float y = n.x;
        float z = n.y;

        // ============================================================
        // Real Spherical Harmonics up to l = 2
        // Convention: Y_l^m, index = l(l + 1) + m
        // Coordinate system: (x, y, z)
        // ============================================================
        //
        //  Index layout table:
        //
        //   l |   m   | index |       Expression
        //  ----|-------|--------|------------------------------------
        //   0 |   0   |   0    | Y_0^0  = 0.282095
        //   1 |  -1   |   1    | Y_1^-1 = 0.488603 * y
        //   1 |   0   |   2    | Y_1^0  = 0.488603 * z
        //   1 |   1   |   3    | Y_1^1  = 0.488603 * x
        //   2 |  -2   |   4    | Y_2^-2 = 1.092548 * xy
        //   2 |  -1   |   5    | Y_2^-1 = 1.092548 * yz
        //   2 |   0   |   6    | Y_2^0  = 0.946175 * z^2 - 0.315392
        //   2 |   1   |   7    | Y_2^1  = 1.092548 * xz
        //   2 |   2   |   8    | Y_2^2  = 0.546274 * (x^2 - y^2)
        // ============================================================

        // with Condon–Shortley phase(-1) ^ m to match Iwasaki sensei code convention
        //  Band 0
        float Y0 = 0.2820947917738781f; // l=0, m=0

        // Band 1
        float Y1 = -0.4886025119029200f * y; // l=1, m=-1
        float Y2 = 0.4886025119029200f * z;  // l=1, m=0
        float Y3 = -0.4886025119029200f * x; // l=1, m=1

        // Band 2
        float Y4 = 1.0925484305920792f * x * y;                       // l=2, m=-2
        float Y5 = -1.0925484305920792f * y * z;                      // l=2, m=-1
        float Y6 = 0.9461746957575601f * z * z - 0.3153915652525200f; // l=2, m=0
        float Y7 = -1.0925484305920792f * x * z;                      // l=2, m=1
        float Y8 = 0.546274f * (x * x - y * y);                       // l=2, m=2

        // Cosine lobe coefficients for irradiance (Peter-Pike Sloan 2002)
        float A[3] = { 3.141593f, 2.094395f, 0.785398f };

        // Dot product: coeffs · basis
        float3 result =
            shCoeffs[0] * Y0 * A[0] +
            shCoeffs[1] * Y1 * A[1] +
            shCoeffs[2] * Y2 * A[1] +
            shCoeffs[3] * Y3 * A[1] +
            shCoeffs[4] * Y4 * A[2] +
            shCoeffs[5] * Y5 * A[2] +
            shCoeffs[6] * Y6 * A[2] +
            shCoeffs[7] * Y7 * A[2] +
            shCoeffs[8] * Y8 * A[2];

        // Clamp negative lighting
        // return max(result, float3(0.f));
        return result;
    }

    static float3 float4ToFloat3_CPU(const float4& v)
    {
        return float3(v.x, v.y, v.z);
    }
}

bool UniformProbeVolume::fetchHermiteSH_CPU(float3 posW, std::vector<float3>& outCoeffs) const
{
    outCoeffs.clear();
    outCoeffs.resize(9, float3(0.0f));

    if (mProbeData.empty())
        return false;

    float3 rel = posW - mMinPoint;

    int3 cell = int3(
        int(std::floor(rel.x / mCellSize.x)),
        int(std::floor(rel.y / mCellSize.y)),
        int(std::floor(rel.z / mCellSize.z))
    );

    int3 maxCell = int3(
        int(mProbeResolution.x) - 1,
        int(mProbeResolution.y) - 1,
        int(mProbeResolution.z) - 1
    );

    cell.x = std::max(0, std::min(cell.x, maxCell.x));
    cell.y = std::max(0, std::min(cell.y, maxCell.y));
    cell.z = std::max(0, std::min(cell.z, maxCell.z));

    float3 cellOrigin = mMinPoint + float3(cell) * mCellSize;

    float3 t = (posW - cellOrigin) / mCellSize;
    t.x = clamp01CPU(t.x);
    t.y = clamp01CPU(t.y);
    t.z = clamp01CPU(t.z);

    auto index = [&](int3 p) -> uint32_t
        {
            return uint32_t(
                p.z * int(mCornerResolution.y) * int(mCornerResolution.x) +
                p.y * int(mCornerResolution.x) +
                p.x
            );
        };

    // Corner order matches your adaptive corner convention:
    // bit 2 = x, bit 1 = y, bit 0 = z
    uint32_t cIdx[8];
    cIdx[0] = index(cell + int3(0, 0, 0));
    cIdx[1] = index(cell + int3(0, 0, 1));
    cIdx[2] = index(cell + int3(0, 1, 0));
    cIdx[3] = index(cell + int3(0, 1, 1));
    cIdx[4] = index(cell + int3(1, 0, 0));
    cIdx[5] = index(cell + int3(1, 0, 1));
    cIdx[6] = index(cell + int3(1, 1, 0));
    cIdx[7] = index(cell + int3(1, 1, 1));

    for (int i = 0; i < 8; ++i)
    {
        if (cIdx[i] >= mProbeData.size())
            return false;
    }

    for (int band = 0; band < 9; ++band)
    {
        float3 v[8];

        // Gradients in Falcor world axes.
        float3 gX[8];
        float3 gY[8];
        float3 gZ[8];

        for (int i = 0; i < 8; ++i)
        {
            const UniformGridCorner& c = mProbeData[cIdx[i]];

            v[i] = float4ToFloat3_CPU(c.coeffs[band]);

            //grad is stored in polar coord
            // Therefore:
            // world X = .y
            // world Y = .z
            // world Z = .x
            gX[i] = float3(c.gradR[band].y, c.gradG[band].y, c.gradB[band].y);
            gY[i] = float3(c.gradR[band].z, c.gradG[band].z, c.gradB[band].z);
            gZ[i] = float3(c.gradR[band].x, c.gradG[band].x, c.gradB[band].x);
        }

        // ------------------------------------------------------------
        // Tricubic Hermite, separable.
        //
        // Corner order:
        // 0=(0,0,0), 1=(0,0,1), 2=(0,1,0), 3=(0,1,1),
        // 4=(1,0,0), 5=(1,0,1), 6=(1,1,0), 7=(1,1,1)
        // ------------------------------------------------------------

        // Pass 1: interpolate along Z.
        float3 q[4];
        float3 q_dX[4];
        float3 q_dY[4];

        for (int i = 0; i < 4; ++i)
        {
            int i0 = 2 * i;
            int i1 = 2 * i + 1;

            q[i] = hermite1D_CPU(
                t.z,
                v[i0],
                gZ[i0] * mCellSize.z,
                v[i1],
                gZ[i1] * mCellSize.z
            );

            // For the next stages, linearly propagate transverse gradients.
            q_dX[i] = lerp(gX[i0], gX[i1], t.z);
            q_dY[i] = lerp(gY[i0], gY[i1], t.z);
        }

        // Pass 2: interpolate along Y.
        float3 r[2];
        float3 r_dX[2];

        for (int i = 0; i < 2; ++i)
        {
            int i0 = 2 * i;
            int i1 = 2 * i + 1;

            r[i] = hermite1D_CPU(
                t.y,
                q[i0],
                q_dY[i0] * mCellSize.y,
                q[i1],
                q_dY[i1] * mCellSize.y
            );

            r_dX[i] = lerp(q_dX[i0], q_dX[i1], t.y);
        }

        // Pass 3: interpolate along X.
        outCoeffs[band] = hermite1D_CPU(
            t.x,
            r[0],
            r_dX[0] * mCellSize.x,
            r[1],
            r_dX[1] * mCellSize.x
        );
    }

    return true;
}

float3 UniformProbeVolume::evaluateIrradianceHermiteCPU(float3 posW, float3 normalW) const
{
    std::vector<float3> coeffs;

    if (!fetchHermiteSH_CPU(posW, coeffs))
        return float3(0.0f);

    return evaluateIrradiance(normalW, coeffs);
}
