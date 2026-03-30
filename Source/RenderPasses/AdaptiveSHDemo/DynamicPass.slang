import Scene.Raster;
import Rendering.Lights.EmissiveLightSampler;
import Rendering.Lights.EmissiveLightSamplerHelpers;
import Utils.Sampling.UniformSampleGenerator;
import Scene.RaytracingInline;
import Utils.Geometry.GeometryHelpers;
#include "Utils/Math/MathConstants.slangh"

// Lightmaps
SamplerState gLinearSampler;
static const bool kUseEmissiveLights = USE_EMISSIVE_LIGHTS;
static const bool kUseAnalyticLights = USE_ANALYTIC_LIGHTS;

cbuffer PerFrameCB
{
    uint gFirstDynamicInstanceID;
    uint sampleIndex;
    EmissiveLightSampler emissiveSampler;
};

VSOut vsMain(VSIn vIn)
{
    VSOut vsOut = defaultVS(vIn);
    return vsOut;
}

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

float3 evalDynamicDirect(ShadingData sd, inout UniformSampleGenerator sg)
{
    //if (!kUseEmissiveLights) return float3(0);

    let lod = ExplicitLodTextureSampler(0.f);
    let mi = gScene.materials.getMaterialInstance(sd, lod);

    PathVertex vertex = PathVertex(0, sd.posW, sd.frame.N, sd.frontFacing);

    int emissiveTriangleCount = gScene.lightCollection.getActiveTriangleCount();
    if (emissiveTriangleCount == 0) return float3(0);

    const int numEmissiveSamplesPerTriangle = 1;
    const int numEmissiveSamples = max(1, numEmissiveSamplesPerTriangle * emissiveTriangleCount);

    float3 direct = float3(0);

    for (int s = 0; s < numEmissiveSamples; ++s)
    {
        LightSample ls;
        if (generateEmissiveSample(vertex, true, sg, ls) && ls.pdf > 0.f)
        {
            SceneRayQuery<false> sceneRayQuery;
            if (sceneRayQuery.traceVisibilityRay(ls.getVisibilityRay()))
            {
                direct += mi.eval(sd, ls.dir, sg)*ls.Li;
            }
        }
    }

    return direct / float(numEmissiveSamples);
}

struct VisibilityRayData
{
    bool visible;
};

float4 psMain(VSOut vsOut, uint triangleIndex : SV_PrimitiveID) : SV_TARGET
{
    uint iid = vsOut.instanceID.index;
   if (iid < gFirstDynamicInstanceID) discard;

    float3 viewDir = normalize(gScene.camera.getPosition() - vsOut.posW);
    ShadingData shadingData = prepareShadingData(vsOut, triangleIndex, viewDir);
    BasicMaterialData baseColorMat = gScene.materials.getBasicMaterialData(vsOut.materialID);

    UniformSampleGenerator sg = UniformSampleGenerator(uint2(vsOut.posH.xy), 0);

    float3  litColor = evalDynamicDirect(shadingData, sg);

    return float4(litColor, 1.0f);
}
