/***************************************************************************
 # LightmapBaking.rt.slang - Ray-Traced Lightmap Baking for Static Objects
 **************************************************************************/

import Scene.Raytracing;
import Utils.Sampling.UniformSampleGenerator;
import Rendering.Lights.LightHelpers;
import Rendering.Lights.EmissiveLightSampler;
import Rendering.Lights.EmissiveLightSamplerHelpers;
import Utils.Geometry.GeometryHelpers;
#include "Utils/Math/MathConstants.slangh"
#include "BakeDataStructures.slang"

// Bindings (adjust registers as per your Falcor pass)
StructuredBuffer<TexelSample> gTexelSamples;
RWStructuredBuffer<float4> gIrradianceAccum;
TextureCube<float3> gEnvMap;
SamplerState gEnvSampler;

static const bool kUseEmissiveLights = USE_EMISSIVE_LIGHTS;
static const bool kUseAnalyticLights = USE_ANALYTIC_LIGHTS;
static const uint kMaxBounces = 3;

cbuffer PerFrameCB
{
    uint sampleIndex;       // Current sample (0 to totalSamples-1)
    uint totalSamples;
    uint numTexels;
    float bias;             // Surface offset bias (e.g., 0.02)
    EmissiveLightSampler emissiveSampler;
};

// Simplified payload (no probe-specific data)
struct PrimaryRayData
{
    float3 Li;
    float3 nextRayOrigin;
    float3 nextRayDir;
    float3 throughput;
    bool terminated;
    UniformSampleGenerator sg;
    uint depth;
};

struct LightSample
{
    float3 Li;
    float pdf;
    float3 origin;
    float distance;
    float3 dir;
    uint lightType;
    Ray getVisibilityRay() { return Ray(origin, dir, 0.f, distance); }
};

struct PathVertex
{
    uint index;
    float3 pos;
    float3 faceNormal;
    bool frontFacing;
    __init(uint index, float3 pos, float3 faceNormal = float3(0.f), bool frontFacing = true)
    {
        this.index = index;
        this.pos = pos;
        this.faceNormal = faceNormal;
        this.frontFacing = frontFacing;
    }
    float3 getRayOrigin(float3 rayDir)
    {
        return computeRayOrigin(pos, dot(faceNormal, rayDir) >= 0 ? faceNormal : -faceNormal);
    }
    float3 getOrientedFaceNormal()
    {
        return frontFacing ? faceNormal : -faceNormal;
    }
};

bool generateEmissiveSample(const PathVertex vertex, const bool upperHemisphere, inout UniformSampleGenerator sg, out LightSample ls)
{
    ls = {};
    if (!kUseEmissiveLights) return false;
    TriangleLightSample tls;
    if (!emissiveSampler.sampleLight(vertex.pos, vertex.getOrientedFaceNormal(), upperHemisphere, sg, tls)) return false;
    ls.Li = tls.pdf > 0.f ? tls.Le / tls.pdf : float3(0);
    ls.pdf = tls.pdf;
    float3 lightPos = computeRayOrigin(tls.posW, tls.normalW);
    ls.origin = vertex.getRayOrigin(lightPos - vertex.pos);
    float3 toLight = lightPos - ls.origin;
    ls.distance = length(toLight);
    ls.dir = normalize(toLight);
    return any(ls.Li > 0.f);
}

// Helpers for cosine-weighted sampling
float3 SampleCosineHemisphere(float2 u)
{
    float phi = 2.0f * M_PI * u.y;
    float cosTheta = sqrt(max(0.0001f, u.x));
    float sinTheta = sqrt(1.0f - u.x);
    return float3(sinTheta * cos(phi), sinTheta * sin(phi), cosTheta);
}

float3 TangentToWorld(float3 localDir, float3 normal)
{
    float3 arb = abs(normal.y) < 0.999f ? float3(0, 1, 0) : float3(1, 0, 0);
    float3 tangent = normalize(cross(arb, normal));
    float3 bitangent = cross(normal, tangent);
    return localDir.x * tangent + localDir.y * bitangent + localDir.z * normal;
}

struct VisibilityRayData
{
    bool visible;
};

// Visibility trace (unchanged)
bool traceVisibilityRay(const Ray ray)
{
    VisibilityRayData hitData;
    hitData.visible = false; // Default to blocked

    RayDesc rayDesc;
    rayDesc.Origin = ray.origin;
    rayDesc.Direction = ray.dir;
    rayDesc.TMin = ray.tMin;
    rayDesc.TMax = ray.tMax;

    TraceRay(
        gScene.rtAccel,
        RAY_FLAG_SKIP_CLOSEST_HIT_SHADER | RAY_FLAG_ACCEPT_FIRST_HIT_AND_END_SEARCH,
        0xFF,
             1, // Ray Index 1 (Shadow/Visibility)
             rayTypeCount,
             1, // Miss Index 1 (Shadow/Visibility)
             rayDesc,
             hitData
    );

    return hitData.visible;
}

[shader("miss")]
void shadowMiss(inout VisibilityRayData hitData)
{
    hitData.visible = true; // Missed all geometry -> Light is visible
}

[shader("anyhit")]
void shadowAnyHit(inout VisibilityRayData hitData, BuiltInTriangleIntersectionAttributes attribs)
{
    // Hit something -> Light is blocked
    hitData.visible = false;
    AcceptHitAndEndSearch();
}

// Miss: Sample env map for background
[shader("miss")]
void primaryMiss(inout PrimaryRayData hitData)
{
    float3 rayDir = normalize(WorldRayDirection());
    // hitData.Li = gEnvMap.SampleLevel(gEnvSampler, rayDir, 0.0f);
    hitData.Li = 0.0f;
    hitData.terminated = true;

    // Optional: * env intensity if in CB
}

// Closest hit: Compute contribution (direct at depth 0, indirect at >0)
[shader("closesthit")]
void primaryClosestHit(inout PrimaryRayData hitData, BuiltInTriangleIntersectionAttributes attribs)
{
    // 1. Geometry & Material Setup
    float3 rayDirW = WorldRayDirection();
    float hitT = RayTCurrent();
    uint triangleIndex = PrimitiveIndex();
    const GeometryInstanceID instanceID = getGeometryInstanceID();
    VertexData v = getVertexData(instanceID, triangleIndex, attribs);
    uint materialID = gScene.getMaterialID(instanceID);
    ShadingData sd = gScene.materials.prepareShadingData(v, materialID, -rayDirW);

    let lod = ExplicitLodTextureSampler(0.f);
    let mi = gScene.materials.getMaterialInstance(sd, lod);
    let bsdfProperties = mi.getProperties(sd);

    // Sample Generator for this hit
    uint3 launchIndex = DispatchRaysIndex();
    // -------------------------------------------------------------------------
    // 2. Direct Lighting (Next Event Estimation / Explicit Light Sampling)
    // -------------------------------------------------------------------------
    // For bounces (depth > 0), we rely on NEE (sampleLight) to handle lighting smoothly.
    float3 accumLight =  bsdfProperties.emission;
   
    // ANALYTIC LIGHTS
    if (kUseAnalyticLights) {
        for (int i = 0; i < gScene.getLightCount(); i++) {
            AnalyticLightSample ls;
            if (sampleLight(sd.posW, gScene.getLight(i), hitData.sg, ls)) {
                // If the light is "above" the face, push out. If "below", push in.
                float3 shadowRayOrigin = computeRayOrigin(sd.posW, dot(sd.faceN, ls.dir) >= 0 ? sd.faceN : -sd.faceN);
                Ray shadowRay;
                shadowRay.origin = shadowRayOrigin; // Offset to avoid acne
                shadowRay.dir = ls.dir;
                shadowRay.tMin = 0.001f;               // Bias
                shadowRay.tMax = ls.distance - 0.001f; // Don't hit the light source itself

                if (traceVisibilityRay(shadowRay))
                {
                    accumLight += mi.eval(sd, ls.dir, hitData.sg) * ls.Li;
                }
            }
        }
    }

    // EMISSIVE LIGHTS
    if (kUseEmissiveLights) {
        const int numEmissiveSamplesPerTriangle = 1;
        int emissiveTriangleCount = gScene.lightCollection.getActiveTriangleCount();
        const int numEmissiveSamples = max(1, (numEmissiveSamplesPerTriangle * emissiveTriangleCount));

        float3 accum = float3(0.0f);
        for (int s = 0; s < numEmissiveSamples; ++s) {
            PathVertex vertex = PathVertex(0, sd.posW, sd.frame.N, sd.frontFacing);
            LightSample ls;
            if (generateEmissiveSample(vertex, true, hitData.sg, ls) && ls.pdf > 0.0f) {
                // FIX: Use traceVisibilityRay with ls.getVisibilityRay()
                if (traceVisibilityRay(ls.getVisibilityRay()))
                {
                    float3 bsdf = mi.eval(sd, ls.dir, hitData.sg); // mi.eval return f(wi,wo)(n dot wi)
                    accum += bsdf * ls.Li;
                }
            }
        }
        accumLight += accum / float(numEmissiveSamples);
    }

    hitData.Li = accumLight;

    // -------------------------------------------------------------------------
    // 3. Indirect Lighting (Generate Next Ray)
    // -------------------------------------------------------------------------
    BSDFSample bsdfSample;

    // mi.sample asks the material: "Where does the ray go next?"
    if (mi.sample(sd, hitData.sg, bsdfSample, true))
    {
        // We use faceN (Geometric Normal) to safely offset the ray
        float3 offsetNormal = dot(sd.faceN, bsdfSample.wo) >= 0 ? sd.faceN : -sd.faceN;
        hitData.nextRayOrigin = computeRayOrigin(sd.posW, offsetNormal);
        hitData.nextRayDir = bsdfSample.wo;
        hitData.throughput = bsdfSample.weight; // This is (F * G * D * Cos) / PDF
        hitData.terminated = false;
    }
    else
    {
        // Ray absorbed (e.g. black diffuse)
        hitData.terminated = true;
        hitData.throughput = float3(0.f);
    }
}

// Raygen: Cosine-weighted sampling from texels
[shader("raygeneration")]
void rayGen()
{
    uint texelIdx = DispatchRaysIndex().x;
    if (texelIdx >= numTexels) return;

    TexelSample ts = gTexelSamples[texelIdx];
    float3 pos = ts.posW;
    float3 normal = normalize(ts.normalW);

    UniformSampleGenerator sg = UniformSampleGenerator(uint2(texelIdx, sampleIndex), sampleIndex);

    float2 xi = sampleNext2D(sg);
    float3 localDir = SampleCosineHemisphere(xi);
    float cosTheta = localDir.z;
    if (cosTheta <= 0.0f) return; // Skip backfacing (rare)

    float3 dir = TangentToWorld(localDir, normal);

    float3 accumColor = float3(0.f);
    float3 throughput = float3(1.f);

    RayDesc ray;
    ray.Origin = pos + normal * bias;
    ray.Direction = dir;
    ray.TMin = 0.001f;
    ray.TMax = 10000.f;

    for (uint depth = 0; depth <= kMaxBounces; ++depth)
    {
        PrimaryRayData hitData;
        hitData.terminated = false;
        hitData.sg = sg;
        hitData.depth = depth;

        TraceRay(gScene.rtAccel, 0, 0xFF, 0, rayTypeCount, 0, ray, hitData);
        sg = hitData.sg;
        accumColor += throughput * hitData.Li;

        if (hitData.terminated)
        {
            break;
        }

        throughput *= hitData.throughput;

        ray.Origin = hitData.nextRayOrigin;
        ray.Direction = hitData.nextRayDir;
    }

    // Accumulate radiance (post-process multiplies by PI / totalSamples)
    gIrradianceAccum[ts.texelIndex] += float4(accumColor, 1.0f);
}
