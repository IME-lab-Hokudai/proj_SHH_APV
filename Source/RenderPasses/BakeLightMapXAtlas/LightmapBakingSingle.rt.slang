// Atlas adapter for Falcor's PathTracer. Keep light transport in the original
// module so baking and camera rendering use the same material/path logic.
import Scene.Raytracing;
import Scene.RaytracingInline;
import RenderPasses.PathTracer.PathTracer;
import Utils.Geometry.GeometryHelpers;
#include "Utils/Math/MathConstants.slangh"
#include "BakeDataStructures.slang"

ParameterBlock<PathTracer> gPathTracer;
StructuredBuffer<TexelSample> gTexelSamples;
RWStructuredBuffer<float4> gIrradianceAccum;

cbuffer PerFrameCB
{
    uint sampleIndex;
    uint numTexels;
    uint atlasWidth;
    uint pageIndex;
};

// Match TracePass.rt.slang: alpha-tested visibility, with dielectric
// transmission handled by BSDF paths rather than transparent shadow shortcuts.
struct BakeVisibilityQuery : IVisibilityQuery
{
    bool traceVisibilityRay(const Ray ray)
    {
        SceneRayQuery<kUseAlphaTest> query;
        return query.traceVisibilityRay(ray);
    }
};

float3 sampleBakeHemisphere(float2 u, float3 normal)
{
    float phi = 2.f * M_PI * u.y;
    float r = sqrt(u.x);
    float3 arb = abs(normal.y) < 0.999f ? float3(0, 1, 0) : float3(1, 0, 0);
    float3 tangent = normalize(cross(arb, normal));
    float3 bitangent = cross(normal, tangent);
    return r * cos(phi) * tangent + r * sin(phi) * bitangent
        + sqrt(max(0.f, 1.f - u.x)) * normal;
}

// The receiving texel is a virtual white Lambertian surface. Its albedo is
// applied by the demo, not baked into the atlas. Store E/pi, as before.
void generateBakePath(const TexelSample texel, out PathState path)
{
    path = {};
    path.setActive();
    path.setVertexIndex(1);
    path.setDiffusePrimaryHit();
    path.setBounces(BounceType::Diffuse, 1);
    path.thp = float3(1.f);
    const uint2 pixel = uint2(texel.texelIndex % atlasWidth, texel.texelIndex / atlasWidth);
    path.id = pixel.x | (pixel.y << 12);
    path.sg = SampleGenerator(pixel, sampleIndex ^ (pageIndex * 0x9e3779b9u));

    const float3 normal = normalize(texel.normalW);
    const float3 faceNormal = normalize(texel.faceNormalW);
    // Falcor uses the oriented geometric normal for light-sampling/MIS state.
    path.normal = faceNormal;
    path.origin = computeRayOrigin(texel.posW, faceNormal);

    // NEE at the receiver is necessary for analytic lights, which cannot be
    // reached by a randomly sampled hemisphere ray. Use Falcor's light sampler
    // and MIS for emissive geometry and the environment as well.
    BakeVisibilityQuery visibility;
    // Supply the actual surface position. Falcor offsets visibility rays itself.
    PathTracer::PathVertex vertex = PathTracer::PathVertex(1, texel.posW, faceNormal, true);
    PathTracer::LightSample ls;
    path.setLightSampled(true, false);
    if (gPathTracer.generateLightSample(vertex, true, false, path.sg, ls))
    {
        const float cosinePdf = max(0.f, dot(normal, ls.dir)) / M_PI;
        float misWeight = 1.f;
        if (ls.lightType != uint(PathTracer::LightType::Analytic))
            misWeight = gPathTracer.evalMIS(1, ls.pdf, 1, cosinePdf);
        if (cosinePdf > 0.f && dot(faceNormal, ls.dir) > 0.f &&
            visibility.traceVisibilityRay(ls.getVisibilityRay()))
            path.L += cosinePdf * misWeight * ls.Li;
    }

    path.dir = sampleBakeHemisphere(sampleNext2D(path.sg), normal);
    path.pdf = max(0.f, dot(normal, path.dir)) / M_PI;
    // Reflection must remain on the receiver's geometric front side. Reject
    // below-surface samples without resampling, preserving the cosine PDF.
    if (path.pdf <= 0.f || dot(faceNormal, path.dir) <= 0.f) path.terminate();
}

[shader("raygeneration")]
void rayGen()
{
    const uint texelIdx = DispatchRaysIndex().x;
    if (texelIdx >= numTexels) return;
    const TexelSample texel = gTexelSamples[texelIdx];

    PathState path;
    generateBakePath(texel, path);
    SceneRayQuery<kUseAlphaTest> sceneQuery;
    BakeVisibilityQuery visibility;

    // Only scheduling differs from Falcor's TracePass. PathState, hit/miss
    // handling, BSDF sampling, nested dielectrics, absorption, NEE, MIS and
    // termination are all supplied by the original PathTracer module.
    while (path.isActive())
    {
        path.incrementVertexIndex();
        float hitT;
        const HitInfo hit = sceneQuery.traceRay(path.getScatterRay(), hitT);
        if (hit.isValid())
        {
            path.setHit(hit);
            gPathTracer.handleHit(path, visibility);
        }
        else
        {
            path.clearHit();
            gPathTracer.handleMiss(path);
        }
    }

    gIrradianceAccum[texel.texelIndex] += float4(path.L, 1.f);
}
