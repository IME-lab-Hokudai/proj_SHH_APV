// Geometry-guided Gaussian blur. Input and output must be different textures.
Texture2D<float4> gInput;
Texture2D<float4> gPosW;
Texture2D<float4> gNormW;
Texture2D<float4> gFilterGuide;
RWTexture2D<float4> gOutput;

cbuffer CB
{
    uint gRadius;
};

[numthreads(16, 16, 1)]
void main(uint3 id : SV_DispatchThreadID)
{
    uint2 res;
    gOutput.GetDimensions(res.x, res.y);
    if (any(id.xy >= res)) return;

    float4 centerPos = gPosW[id.xy];
    if (centerPos.w <= 0.f)
    {
        gOutput[id.xy] = float4(0.f);
        return;
    }

    float3 center = gInput[id.xy].rgb;
    if (gRadius == 0)
    {
        gOutput[id.xy] = float4(center, 1.f);
        return;
    }

    float3 normal = normalize(gNormW[id.xy].xyz);
    float4 guide = gFilterGuide[id.xy];
    int radius = int(min(gRadius, 8u));
    float sigma = max(0.5f * float(radius), 0.5f);
    float3 sum = center;
    float weightSum = 1.f;

    for (int y = -radius; y <= radius; ++y)
    {
        for (int x = -radius; x <= radius; ++x)
        {
            if (x == 0 && y == 0) continue;
            int2 p = int2(id.xy) + int2(x, y);
            if (any(p < int2(0)) || any(p >= int2(res))) continue;

            float4 pos = gPosW[p];
            if (pos.w <= 0.f) continue;
            if (any(gFilterGuide[p].zw != guide.zw)) continue;

            float alignment = dot(normal, normalize(gNormW[p].xyz));
            if (alignment < 0.9f) continue;

            // Reject UV neighbors far apart in the scene. The threshold scales
            // with this triangle's texel footprint, not scene units or orientation.
            float reach = dot(abs(float2(x, y)), guide.xy);
            float3 delta = pos.xyz - centerPos.xyz;
            float maxDistance = max(1.5f * reach, 1e-6f);
            if (dot(delta, delta) > maxDistance * maxDistance) continue;

            float weight = exp(-float(x * x + y * y) / (2.f * sigma * sigma));
            weight *= pow(saturate(alignment), 16.f);
            sum += weight * gInput[p].rgb;
            weightSum += weight;
        }
    }

    // Do not reject neighbors by brightness: real black samples remain valid.
    gOutput[id.xy] = float4(sum / weightSum, 1.f);
}
