/***************************************************************************
 # Copyright (c) 2015-24, NVIDIA CORPORATION. All rights reserved.
 #
 # Redistribution and use in source and binary forms, with or without
 # modification, are permitted provided that the following conditions
 # are met:
 #  * Redistributions of source code must retain the above copyright
 #    notice, this list of conditions and the following disclaimer.
 #  * Redistributions in binary form must reproduce the above copyright
 #    notice, this list of conditions and the following disclaimer in the
 #    documentation and/or other materials provided with the distribution.
 #  * Neither the name of NVIDIA CORPORATION nor the names of its
 #    contributors may be used to endorse or promote products derived
 #    from this software without specific prior written permission.
 #
 # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS "AS IS" AND ANY
 # EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 # IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 # PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 # CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 # EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 # PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 # PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 # OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 # (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 # OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 **************************************************************************/
import Scene.Raytracing;
import Utils.Sampling.UniformSampleGenerator;
import Rendering.Lights.LightHelpers;
import ProbeSamplingData;
import Rendering.Lights.EmissiveLightSampler;
import Rendering.Lights.EmissiveLightSamplerHelpers;
import Utils.Geometry.GeometryHelpers;

StructuredBuffer<ProbeDirSample> gProbeDirSamples;
StructuredBuffer<float3> gProbePositions;
RWStructuredBuffer<float4> gProbeSamplingOutput;
EmissiveLightSampler emissiveSampler;

static const bool kUseEmissiveLights = USE_EMISSIVE_LIGHTS;
static const bool kUseAnalyticLights = USE_ANALYTIC_LIGHTS;

cbuffer PerFrameCB
{
    uint sampleIndex; // for UniformSampleGenerator. I don't know how to use it properly yet'
    uint numSamplePerProbe;
};

struct PrimaryRayData
{
    float4 color;
    uint depth;
    float hitT;
};

/** Describes a light sample.
Copy from Falcor's PathTracer.slang
 */
struct LightSample
{
    float3 Li;      ///< Incident radiance at the shading point (unshadowed). This is already divided by the pdf.
    float pdf;      ///< Pdf with respect to solid angle at the shading point.
    float3 origin;  ///< Ray origin for visibility evaluation (offseted to avoid self-intersection).
    float distance; ///< Ray distance for visibility evaluation (shortened to avoid self-intersection).
    float3 dir;     ///< Ray direction for visibility evaluation (normalized).
    uint lightType; ///< Light type this sample comes from (LightType casted to uint).

    Ray getVisibilityRay() { return Ray(origin, dir, 0.f, distance); }
};

/** Describes a path vertex.
Copy from Falcor's PathTracer.slang
 */
struct PathVertex
{
    uint index;        ///< Vertex index (0 = camera, 1 = primary hit, 2 = secondary hit, etc.).
    float3 pos;        ///< Vertex position.
    float3 faceNormal; ///< Geometry normal at the vertex (zero if not on a surface).
    bool frontFacing;  ///< True if path vertex is on the front-facing side (if on a surface).

    /** Initializes a path vertex.
        \param[in] index Vertex index.
        \param[in] pos Vertex position.
        \param[in] faceNormal Geometry normal.
        \param[in] frontFacing Front-facing flag.
    */
    __init(uint index, float3 pos, float3 faceNormal = float3(0.f), bool frontFacing = true)
    {
        this.index = index;
        this.pos = pos;
        this.faceNormal = faceNormal;
        this.frontFacing = frontFacing;
    }

    /** Get position with offset applied in direction of the geometry normal to avoid self-intersection
        for visibility rays.
        \param[in] rayDir Direction of the visibility ray (does not need to be normalized).
        \return Returns the offseted position.
    */
    float3 getRayOrigin(float3 rayDir)
    {
        return computeRayOrigin(pos, dot(faceNormal, rayDir) >= 0 ? faceNormal : -faceNormal);
    }

    /** Returns the oriented face normal.
        \return Face normal flipped to the same side as the view vector.
    */
    float3 getOrientedFaceNormal()
    {
        return frontFacing ? faceNormal : -faceNormal;
    }
};

/** Generates a light sample on the emissive geometry.
     \param[in] vertex Path vertex.
     \param[in] upperHemisphere True if only upper hemisphere should be considered.
     \param[in,out] sg Sample generator.
     \param[out] ls Struct describing valid samples.
     \return True if the sample is valid and has nonzero contribution, false otherwise.
 */
bool generateEmissiveSample(const PathVertex vertex, const bool upperHemisphere, inout UniformSampleGenerator sg, out LightSample ls)
{
    ls = {}; // Default initialization to avoid divergence at returns.
    if (!kUseEmissiveLights) return false;

    TriangleLightSample tls;
    if (!emissiveSampler.sampleLight(vertex.pos, vertex.getOrientedFaceNormal(), upperHemisphere, sg, tls)) return false;

    // Setup returned sample.
    ls.Li = tls.pdf > 0.f ? tls.Le / tls.pdf : float3(0);
    ls.pdf = tls.pdf;
    // Offset shading and light position to avoid self-intersection.
    float3 lightPos = computeRayOrigin(tls.posW, tls.normalW);
    ls.origin = vertex.getRayOrigin(lightPos - vertex.pos);
    float3 toLight = lightPos - ls.origin;
    ls.distance = length(toLight);
    ls.dir = normalize(toLight);

    return any(ls.Li > 0.f);
}

//struct ShadowRayData
//{
//    bool hit;
//};

//[shader("miss")]
//void shadowMiss(inout ShadowRayData hitData)
//{
//    hitData.hit = false;
//}

//[shader("anyhit")]
//void shadowAnyHit(inout ShadowRayData hitData, BuiltInTriangleIntersectionAttributes attribs)
//{
//    hitData.hit = true;
//}

[shader("miss")]
void primaryMiss(inout PrimaryRayData hitData)
{
    float3 rayOrgW = WorldRayOrigin();

    float3 color = float3(0.0f, 0.0f, 0.0f);
    // Direct lighting from analytic light sources.
    for (int i = 0; i < gScene.getLightCount(); i++)
    {
        AnalyticLightSample ls;
        if (evalLightApproximate(rayOrgW, gScene.getLight(i), ls)) // there is another function that handle more type of light
        {
            color += ls.Li; // when miss we can take light color directly
        }
    }
    hitData.color = float4(color, 1.0f);
    hitData.hitT = -1.f;
}

//bool traceShadowRay(RayDesc ray)
//{
//    ShadowRayData rayData;
//    rayData.hit = true;
//    TraceRay(gScene.rtAccel, RAY_FLAG_ACCEPT_FIRST_HIT_AND_END_SEARCH, 0xFF, 1 /* ray index 1 is shadow ray */, rayTypeCount, 1, ray, rayData);
//    return rayData.hit;
//}

//float3 getReflectionColor(const ShadingData sd, const BSDFProperties bsdfProperties, float3 normalW, float3 worldRayDir, uint hitDepth)
//{
//    float3 reflectColor = float3(0, 0, 0);
//    if (hitDepth == 0)
//    {
//        PrimaryRayData secondaryRay = {};
//        secondaryRay.depth = 1;

//        RayDesc ray;
//        ray.Origin = sd.posW;
//        ray.Direction = reflect(worldRayDir, normalW);
//        ray.TMin = 0.001f;
//        ray.TMax = 100000.f;
//        TraceRay(gScene.rtAccel, 0 /*rayFlags*/, 0xFF, 0 /* ray index*/, rayTypeCount, 0, ray, secondaryRay);

//        // A very non-PBR inaccurate way to do reflections.
//        reflectColor = secondaryRay.hitT < 0.f ? float3(0.f) : secondaryRay.color.rgb;
//        float falloff = max(1.f, (secondaryRay.hitT * secondaryRay.hitT));
//        reflectColor *= 20.f / falloff;
//        float ggxAlpha = min(0.5f, max(1e-8f, bsdfProperties.roughness * bsdfProperties.roughness));
//        reflectColor += bsdfProperties.specularReflectance * reflectColor * (ggxAlpha * ggxAlpha);
//    }
//    return reflectColor;
//}

[shader("closesthit")]
void primaryClosestHit(inout PrimaryRayData hitData, BuiltInTriangleIntersectionAttributes attribs)
{
    // Get the hit-point data.
    float3 rayDirW = WorldRayDirection();
    float hitT = RayTCurrent();
    uint triangleIndex = PrimitiveIndex();

    // Prepare the shading data.
    const GeometryInstanceID instanceID = getGeometryInstanceID();
    VertexData v = getVertexData(instanceID, triangleIndex, attribs);
    uint materialID = gScene.getMaterialID(instanceID);
    ShadingData sd = gScene.materials.prepareShadingData(v, materialID, -rayDirW);

    // Create material instance and query its properties.
    let lod = ExplicitLodTextureSampler(0.f);
    let mi = gScene.materials.getMaterialInstance(sd, lod); // IMaterialInstance
    let bsdfProperties = mi.getProperties(sd);

    // Get emission (self-emissive surfaces)
     float3 color = bsdfProperties.emission;

     uint3 launchIndex = DispatchRaysIndex();

     UniformSampleGenerator sg = UniformSampleGenerator(launchIndex.xy, sampleIndex);

     // Direct lighting from analytic light sources.
    if (kUseAnalyticLights) {
         for (int i = 0; i < gScene.getLightCount(); i++)
         {
              AnalyticLightSample ls;
              if (evalLightApproximate(sd.posW, gScene.getLight(i), ls)) // there is another function that handle more type of light
              {
                  color += mi.eval(sd, ls.dir, sg) * ls.Li;
              }
         }
    }

    if(kUseEmissiveLights){
        const int numEmissiveSamplesPerTriangle = 4; // hard coded for now
        int emissiveTriangleCount = gScene.lightCollection.getActiveTriangleCount();
        const int numEmissiveSamples = max(1, (numEmissiveSamplesPerTriangle * emissiveTriangleCount)); // at least 1 sample per triangle
        float3 accum = float3(0.0f);
        for (int s = 0; s < numEmissiveSamples; ++s)
        {
            // Prepare a PathVertex for the current shading point
            PathVertex vertex = PathVertex(0, sd.posW, sd.faceN, sd.frontFacing);

            LightSample ls;
            // Sample the upper hemisphere (typical for diffuse/irradiance)
            if (generateEmissiveSample(vertex, true, sg, ls) && ls.pdf > 0.0f)
            {
                float3 bsdf = mi.eval(sd, ls.dir, sg);
                accum += bsdf * ls.Li / ls.pdf;
            }
        }
        color += accum / float(numEmissiveSamples);
    }
   
    hitData.color = float4(color, 1.0f);
    hitData.hitT = hitT;
}

//[shader("anyhit")]
//void primaryAnyHit(inout PrimaryRayData rayData, BuiltInTriangleIntersectionAttributes attribs)
//{
//    // Alpha test for non-opaque geometry.
//    GeometryInstanceID instanceID = getGeometryInstanceID();
//    VertexData v = getVertexData(instanceID, PrimitiveIndex(), attribs);
//    uint materialID = gScene.getMaterialID(instanceID);
//    if (gScene.materials.alphaTest(v, materialID, 0.f))
//    {
//        IgnoreHit();
//    }
//}

[shader("raygeneration")]
void rayGen()
{
    uint dirIdx = DispatchRaysIndex().x;
    uint probeIdx = DispatchRaysIndex().y;

    float3 probePos = gProbePositions[probeIdx];
    float3 dir = normalize(gProbeDirSamples[dirIdx].dir);

    RayDesc ray;
    ray.Origin = probePos;
    ray.Direction = dir;
    ray.TMin = 0.001f;
    ray.TMax = 100.f;

    PrimaryRayData hitData;
    hitData.depth = 0;

    TraceRay(gScene.rtAccel, 0, 0xFF, 0, rayTypeCount, 0, ray, hitData);

    // Store result in a 2D (flattened) buffer
    gProbeSamplingOutput[probeIdx * numSamplePerProbe + dirIdx] = hitData.color;
}
