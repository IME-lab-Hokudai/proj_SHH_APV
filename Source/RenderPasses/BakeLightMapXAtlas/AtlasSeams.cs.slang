import Scene.Raster;
#include "BakeDataStructures.slang"

cbuffer CB
{
    uint gCount;
    uint gVertexOffset;
    uint gRadius;
};

RWStructuredBuffer<float3> gNormals;
StructuredBuffer<uint2> gPixels;
RWStructuredBuffer<SeamTexel> gReadback;
StructuredBuffer<float4> gValues;
Texture2D<float4> gInput;
Texture2D<float4> gNormW;
Texture2D<float4> gPosW;
Texture2D<float4> gFilterGuide;
Texture2D<uint> gTriangleID;
Texture2D<float4> gSeeds;
StructuredBuffer<uint2> gChartRanges;
StructuredBuffer<uint> gChartIDs;
RWTexture2D<float4> gOutput;

[numthreads(64, 1, 1)]
void readNormals(uint3 id : SV_DispatchThreadID)
{
    if (id.x >= gCount) return;
    gNormals[id.x] = gScene.getVertex(gVertexOffset + id.x).normal;
}

[numthreads(64, 1, 1)]
void gather(uint3 id : SV_DispatchThreadID)
{
    if (id.x >= gCount) return;
    uint2 p = gPixels[id.x];
    SeamTexel texel;
    texel.lighting = gInput[p];
    texel.normalCoverage = gNormW[p];
    uint2 instanceParts = uint2(gFilterGuide[p].zw);
    texel.owner = uint2(instanceParts.x | (instanceParts.y << 16), gTriangleID[p]);
    gReadback[id.x] = texel;
}

uint getRegion(uint2 p)
{
    uint2 parts = uint2(gFilterGuide[p].zw);
    uint instance = parts.x | (parts.y << 16);
    uint count, stride;
    gChartRanges.GetDimensions(count, stride);
    if (instance >= count) return 0xffffffffu;
    uint2 range = gChartRanges[instance];
    uint triangle = gTriangleID[p];
    return triangle < range.y ? gChartIDs[range.x + triangle] : 0xffffffffu;
}

[numthreads(16, 16, 1)]
void spreadCorrection(uint3 id : SV_DispatchThreadID)
{
    uint width, height;
    gOutput.GetDimensions(width, height);
    if (id.x >= width || id.y >= height) return;
    uint2 p = id.xy;
    float4 seed = gSeeds[p];
    float4 normal = gNormW[p];
    if (seed.w > 0.f)
    {
        gOutput[p] = float4(max(gInput[p].rgb + seed.rgb, 0.f), 1.f);
        return;
    }
    if (normal.w <= 0.f) { gOutput[p] = float4(0.f); return; }
    uint region = getRegion(p);
    float3 delta = float3(0.f);
    float weightSum = 0.f;
    float nearest = float(gRadius + 1);
    for (int y = -int(gRadius); y <= int(gRadius); ++y)
    for (int x = -int(gRadius); x <= int(gRadius); ++x)
    {
        int2 q = int2(p) + int2(x, y);
        if (any(q < 0) || q.x >= int(width) || q.y >= int(height)) continue;
        float distance = length(float2(x, y));
        if (distance > float(gRadius)) continue;
        float4 neighbor = gSeeds[q];
        if (neighbor.w <= 0.f || gNormW[q].w <= 0.f || region == 0xffffffffu || getRegion(uint2(q)) != region) continue;
        float reach = max(2.f * dot(abs(float2(x, y)), gFilterGuide[p].xy), 1e-6f);
        float3 worldDelta = gPosW[q].xyz - gPosW[p].xyz;
        if (dot(worldDelta, worldDelta) > reach * reach) continue;
        float weight = exp(-distance * distance / max(0.5f * float(gRadius * gRadius), 1.f));
        delta += weight * neighbor.rgb;
        weightSum += weight;
        nearest = min(nearest, distance);
    }
    // Fade the correction, not the original lighting, into the patch interior.
    float fade = saturate(1.f - nearest / float(gRadius + 1));
    fade = fade * fade * (3.f - 2.f * fade);
    float3 color = gInput[p].rgb + fade * delta / max(weightSum, 1e-8f);
    gOutput[p] = float4(max(color, 0.f), 1.f);
}

[numthreads(64, 1, 1)]
void apply(uint3 id : SV_DispatchThreadID)
{
    if (id.x >= gCount) return;
    // CPU deduplicates pixel addresses: each invocation has exclusive ownership.
    gOutput[gPixels[id.x]] = gValues[id.x];
}
