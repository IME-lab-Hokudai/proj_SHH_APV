// Uniform SH lookup for the shared xatlas scene renderer.
// Buffer order is X-fastest. Interpolation corner bits are Z, Y, X.
// Stored gradients use polar axes: world X/Y/Z map to gradient y/z/x.
StructuredBuffer<UniformGridCorner> gUniformProbes;

cbuffer PerFrameCB
{
    float3 gGridMin;
    float3 gCellSize;
    uint3 gProbeCountDim;
};

uint getUniformGridIndex(int3 coord)
{
    uint strideY = gProbeCountDim.x;
    uint strideZ = gProbeCountDim.x * gProbeCountDim.y;

    return coord.x + (coord.y * strideY) + (coord.z * strideZ);
}

void fetchTrilinearInterpolatedSH(float3 pos, out float3 finalCoeffs[9])
{
    float3 relPos = pos - gGridMin;

    int3 cellCoords = int3(floor(relPos / gCellSize));

    int3 maxCellIndex = int3(gProbeCountDim) - int3(1, 1, 1);
    cellCoords = max(int3(0), min(cellCoords, maxCellIndex - int3(1)));

    float3 cellOrigin = gGridMin + float3(cellCoords) * gCellSize;
    float3 t = saturate((pos - cellOrigin) / gCellSize);
    int cIdx[8];
    cIdx[0] = getUniformGridIndex(cellCoords + int3(0, 0, 0)); // 000
    cIdx[1] = getUniformGridIndex(cellCoords + int3(0, 0, 1)); // 001
    cIdx[2] = getUniformGridIndex(cellCoords + int3(0, 1, 0)); // 010
    cIdx[3] = getUniformGridIndex(cellCoords + int3(0, 1, 1)); // 011
    cIdx[4] = getUniformGridIndex(cellCoords + int3(1, 0, 0)); // 100
    cIdx[5] = getUniformGridIndex(cellCoords + int3(1, 0, 1)); // 101
    cIdx[6] = getUniformGridIndex(cellCoords + int3(1, 1, 0)); // 110
    cIdx[7] = getUniformGridIndex(cellCoords + int3(1, 1, 1)); // 111

    [unroll]
    for (int band = 0; band < 9; ++band)
    {
        float3 v[8];
        float3 g_FalcorX[8];
        float3 g_FalcorY[8];
        float3 g_FalcorZ[8];
        [unroll]
        for (int i = 0; i < 8; ++i)
        {
            v[i] = gUniformProbes[cIdx[i]].coeffs[band].rgb;
        }

        float3 z00 = lerp(v[0], v[1], t.z);
        float3 z01 = lerp(v[2], v[3], t.z);
        float3 z10 = lerp(v[4], v[5], t.z);
        float3 z11 = lerp(v[6], v[7], t.z);

        float3 y0 = lerp(z00, z01, t.y);
        float3 y1 = lerp(z10, z11, t.y);

        finalCoeffs[band] = lerp(y0, y1, t.x);
    }
}

void fetchHermiteInterpolatedSH(float3 pos, out float3 finalCoeffs[9])
{

    float3 relPos = pos - gGridMin;

    int3 cellCount = int3(gProbeCountDim) - int3(1, 1, 1);

    int3 cellCoords = int3(floor(relPos / gCellSize));
    cellCoords = clamp(cellCoords, int3(0), cellCount - int3(1, 1, 1));

    float3 cellOrigin = gGridMin + float3(cellCoords) * gCellSize;
    float3 t = (pos - cellOrigin) / gCellSize;
    t = clamp(t, 0.0f, 1.0f);
    int cIdx[8];
    cIdx[0] = getUniformGridIndex(cellCoords + int3(0, 0, 0)); // 000
    cIdx[1] = getUniformGridIndex(cellCoords + int3(0, 0, 1)); // 001
    cIdx[2] = getUniformGridIndex(cellCoords + int3(0, 1, 0)); // 010
    cIdx[3] = getUniformGridIndex(cellCoords + int3(0, 1, 1)); // 011
    cIdx[4] = getUniformGridIndex(cellCoords + int3(1, 0, 0)); // 100
    cIdx[5] = getUniformGridIndex(cellCoords + int3(1, 0, 1)); // 101
    cIdx[6] = getUniformGridIndex(cellCoords + int3(1, 1, 0)); // 110
    cIdx[7] = getUniformGridIndex(cellCoords + int3(1, 1, 1)); // 111

    float3 size = gCellSize;
    for (int band = 0; band < 9; ++band)
    {
        float3 v[8];
        float3 g_FalcorX[8];
        float3 g_FalcorY[8];
        float3 g_FalcorZ[8];

        for (int i = 0; i < 8; ++i)
        {
            UniformGridCorner c = gUniformProbes[cIdx[i]];
            v[i] = c.coeffs[band].rgb;

            g_FalcorX[i] = float3(c.gradR[band].y, c.gradG[band].y, c.gradB[band].y);
            g_FalcorY[i] = float3(c.gradR[band].z, c.gradG[band].z, c.gradB[band].z);
            g_FalcorZ[i] = float3(c.gradR[band].x, c.gradG[band].x, c.gradB[band].x);
        }

        float3 q[4];
        float3 q_dX[4];
        float3 q_dY[4];

        for (int i = 0; i < 4; ++i)
        {
            int i0 = 2 * i;
            int i1 = 2 * i + 1;

            q[i] = hermite1D(t.z, v[i0], g_FalcorZ[i0] * size.z, v[i1], g_FalcorZ[i1] * size.z);

            q_dX[i] = lerp(g_FalcorX[i0], g_FalcorX[i1], t.z);
            q_dY[i] = lerp(g_FalcorY[i0], g_FalcorY[i1], t.z);
        }

        float3 r[2];
        float3 r_dX[2];

        for (int i = 0; i < 2; ++i)
        {
            int i0 = 2 * i;
            int i1 = 2 * i + 1;

            r[i] = hermite1D(t.y, q[i0], q_dY[i0] * size.y, q[i1], q_dY[i1] * size.y);

            r_dX[i] = lerp(q_dX[i0], q_dX[i1], t.y);
        }

        finalCoeffs[band] = hermite1D(t.x, r[0], r_dX[0] * size.x, r[1], r_dX[1] * size.x);
    }
}

bool sampleGridSH(float3 pos, out float3 coeffs[MAX_SH_BASIS])
{
    [unroll]
    for (int band = 0; band < MAX_SH_BASIS; ++band) coeffs[band] = float3(0.f);
    if (any(gProbeCountDim < uint3(2)) || any(gCellSize <= float3(0.f)))
        return false;

    float3 gridMax = gGridMin + float3(gProbeCountDim - uint3(1)) * gCellSize;
    if (any(pos < gGridMin) || any(pos > gridMax))
        return false;

    fetchHermiteInterpolatedSH(pos, coeffs);
    return true;
}
