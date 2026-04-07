Texture2D<float4> gInput;
Texture2D<float>  gMask;
RWTexture2D<float4> gOutput;

cbuffer CB
{
    uint2 gFrameDim;
};

[numthreads(8, 8, 1)]
void main(uint3 tid : SV_DispatchThreadID)
{
    uint2 p = tid.xy;
    if (any(p >= gFrameDim)) return;

    if (gMask[p] < 0.5f)
    {
        gOutput[p] = float4(0,0,0,0);
        return;
    }

    float3 sum = 0;
    float wsum = 0;

    for (int dy = -2; dy <= 2; ++dy)
    {
        for (int dx = -2; dx <= 2; ++dx)
        {
            int2 q = int2(p) + int2(dx, dy);
            if (q.x < 0 || q.y < 0 || q.x >= int(gFrameDim.x) || q.y >= int(gFrameDim.y)) continue;
            if (gMask[q] < 0.5f) continue;

            float w = 1.0f; // start simple
            sum += gInput[q].rgb * w;
            wsum += w;
        }
    }

    gOutput[p] = float4(wsum > 0 ? sum / wsum : gInput[p].rgb, 1.0f);
}
