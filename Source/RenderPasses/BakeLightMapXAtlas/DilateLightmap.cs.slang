// Fill atlas gutters AFTER blur, using only original covered texels as seeds.
Texture2D<float4> gInput;
RWTexture2D<float4> gOutput;

cbuffer CB
{
    uint gPadding;
};

[numthreads(16, 16, 1)]
void main(uint3 id : SV_DispatchThreadID)
{
    uint2 res;
    gOutput.GetDimensions(res.x, res.y);
    if (any(id.xy >= res)) return;

    float4 value = gInput[id.xy];
    if (value.a == 0.f)
    {
        int radius = int(gPadding);
        int bestDistanceSquared = radius * radius + 1;
        for (int y = -radius; y <= radius; ++y)
        {
            for (int x = -radius; x <= radius; ++x)
            {
                int distanceSquared = x * x + y * y;
                if (distanceSquared >= bestDistanceSquared) continue;
                int2 p = int2(id.xy) + int2(x, y);
                if (any(p < int2(0)) || any(p >= int2(res))) continue;
                float4 candidate = gInput[p];
                if (candidate.a > 0.f)
                {
                    value = candidate;
                    bestDistanceSquared = distanceSquared;
                }
            }
        }
    }
    gOutput[id.xy] = float4(value.rgb, 1.f);
}
