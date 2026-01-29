/***************************************************************************
 # ProbeSampling.rt.slang - Updated for 5-Bounce Path Tracing
 **************************************************************************/

import Scene.Raytracing;
import Utils.Sampling.UniformSampleGenerator;
import Rendering.Lights.LightHelpers;
import ProbeSamplingData;
import Rendering.Lights.EmissiveLightSampler;
import Rendering.Lights.EmissiveLightSamplerHelpers;
import Utils.Geometry.GeometryHelpers;
#include "Utils/Math/MathConstants.slangh"

StructuredBuffer<ProbeDirSample> gProbeDirSamples;
StructuredBuffer<float3> gProbePositions;
RWStructuredBuffer<ProbeSampleData> gProbeSamplingOutput;

static const bool kUseEmissiveLights = USE_EMISSIVE_LIGHTS;
static const bool kUseAnalyticLights = USE_ANALYTIC_LIGHTS;
static const uint kMaxBounces = 5; 

cbuffer PerFrameCB
{
    uint sampleIndex;
    uint numSamplePerProbe;
    EmissiveLightSampler emissiveSampler;
};

// Updated Payload to handle "Next Ray" data
struct PrimaryRayData
{
    float3 Li;           // Direct Light contribution at this hit
    float3 nextRayOrigin;// Where the next bounce starts
    float3 nextRayDir;   // Where the next bounce goes
    float3 throughput;   // Color weight (BRDF * Cos / PDF)
    bool terminated;     // Did the ray hit the sky or get absorbed?
    
    // Probe Data (Only saved on 1st hit)
    float3 s;            // Hit Position
    float3 n;            // Hit Normal
    float hitT;

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

struct VisibilityRayData
{
    bool visible;
};

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

[shader("miss")]
void primaryMiss(inout PrimaryRayData hitData)
{
    // Background / Environment Map support could go here.
    // For now, space is black.
    hitData.Li = float3(0.f); 
    hitData.terminated = true; // Stop bouncing
    
    // Safety defaults
    hitData.s = float3(0.f);
    hitData.n = float3(0.f);
    hitData.hitT = -1.f;
}

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
    
    // Save Probe Data (Position/Normal)
    hitData.s = sd.posW;
    hitData.n = normalize(sd.frame.N);
    hitData.hitT = hitT;

    let lod = ExplicitLodTextureSampler(0.f);
    let mi = gScene.materials.getMaterialInstance(sd, lod);
    let bsdfProperties = mi.getProperties(sd);
    
    // Sample Generator for this hit
    uint3 launchIndex = DispatchRaysIndex();
    // We mix hitT into the seed so every bounce gets different random numbers
    //UniformSampleGenerator sg = UniformSampleGenerator(launchIndex.xy, sampleIndex);

    // -------------------------------------------------------------------------
    // 2. Direct Lighting (Next Event Estimation / Explicit Light Sampling)
    // -------------------------------------------------------------------------
    // For bounces (depth > 0), we rely on NEE (sampleLight) to handle lighting smoothly.
     float3 directLight = (hitData.depth == 0) ? bsdfProperties.emission : float3(0.0f);
    //float3 directLight = float3(0.0f);

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
                    directLight += mi.eval(sd, ls.dir, hitData.sg) * ls.Li;
                }
            }
        }
    }

    // EMISSIVE LIGHTS
    if (kUseEmissiveLights) {
        const int numEmissiveSamplesPerTriangle = 4;
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
                    float3 bsdf = mi.eval(sd, ls.dir, hitData.sg);
                    accum += bsdf * ls.Li;
                }
            }
        }
        directLight += accum / float(numEmissiveSamples);
    }
    
    hitData.Li = directLight;

    // -------------------------------------------------------------------------
    // 3. Indirect Lighting (Generate Next Ray)
    // -------------------------------------------------------------------------
    BSDFSample bsdfSample;

    // mi.sample asks the material: "Where does the ray go next?"
    if (mi.sample(sd, hitData.sg, bsdfSample, true))
    {
        // We use faceN (Geometric Normal) to safely offset the ray
       // float3 offsetNormal = dot(sd.faceN, bsdfSample.wo) >= 0 ? sd.faceN : -sd.faceN;
        hitData.nextRayOrigin = computeRayOrigin(sd.posW, sd.faceN);
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

[shader("raygeneration")]
void rayGen()
{
    uint dirIdx = DispatchRaysIndex().x;
    uint probeIdx = DispatchRaysIndex().y;
    
    // 1. Initial Setup
    float3 probePos = gProbePositions[probeIdx];
    float3 dir = gProbeDirSamples[dirIdx].dir;

    // Convert Polar Z-up to Falcor Y-up
    float3 dirFalcor;
    dirFalcor.x = dir.y; 
    dirFalcor.y = dir.z; 
    dirFalcor.z = dir.x; 

    // 2. Path Tracing State
    float3 accumColor = float3(0.f);
    float3 throughput = float3(1.f);
   
    RayDesc ray;
    ray.Origin = probePos;
    ray.Direction = dirFalcor;
    ray.TMin = 0.001f;
    ray.TMax = 1000.f; // Increased range

    // Variables to store the FIRST hit data (for the probe)
    float3 firstHitPos = float3(0.f);
    float3 firstHitNorm = float3(0.f);
    float firstHitT = -1.0f;

    //UniformSampleGenerator sg = UniformSampleGenerator(DispatchRaysIndex().xy, sampleIndex);
    UniformSampleGenerator sg = UniformSampleGenerator(uint2(dirIdx, 0), sampleIndex);
    // 3. The Bounce Loop
    for (uint depth = 0; depth <= kMaxBounces; ++depth)
    {
        PrimaryRayData hitData;
        hitData.terminated = false;
        hitData.sg = sg;
        hitData.depth = depth;
        // Trace
        TraceRay(gScene.rtAccel, 0, 0xFF, 0, rayTypeCount, 0, ray, hitData);
        sg = hitData.sg; // Update sample generator state
        // Capture First Hit Data Only (The probe only cares about the first wall)
        if (depth == 0)
        {
            firstHitPos = hitData.s;
            firstHitNorm = hitData.n;
            firstHitT = hitData.hitT;
        }

        if (hitData.terminated)
        {
            // Add environment light if you had one here
            break; 
        }

        // Accumulate Light: Light = PreviousThroughput * CurrentDirectLight
        accumColor += throughput * hitData.Li;

        // Prepare for next bounce
        throughput *= hitData.throughput;
        
        // Russian Roulette (Optional Optimization)
        // Stops rays that become too dark to matter
        //if (depth > 2)
        //{
        //    float p = max(throughput.r, max(throughput.g, throughput.b));
        //    if (p < 0.01f) break; // Kill ray
        //}

        // Update Ray
        ray.Origin = hitData.nextRayOrigin;
        ray.Direction = hitData.nextRayDir;
    }

    // 4. Store Result (Convert back to Polar coords for storage)
    // Normal
    float3 nFalcor = firstHitNorm;
    float3 nPolar;
    nPolar.x = nFalcor.z; 
    nPolar.y = nFalcor.x; 
    nPolar.z = nFalcor.y; 

    // Position
    float3 sFalcor = firstHitPos;
    float3 sPolar;
    sPolar.x = sFalcor.z;
    sPolar.y = sFalcor.x; 
    sPolar.z = sFalcor.y; 

    uint outIdx = probeIdx * numSamplePerProbe + dirIdx;
    gProbeSamplingOutput[outIdx].Li = float4(accumColor, 1.0f); // Final Indirect Color
    gProbeSamplingOutput[outIdx].n = float4(nPolar, 1.0f);
    gProbeSamplingOutput[outIdx].s = float4(sPolar, 1.0f);
    gProbeSamplingOutput[outIdx].hitT = firstHitT;
}
