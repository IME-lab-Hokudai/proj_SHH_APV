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

//float3 evalDynamicDirect(ShadingData sd, inout UniformSampleGenerator sg)
//{
//    //if (!kUseEmissiveLights) return float3(0);

//    let lod = ExplicitLodTextureSampler(0.f);
//    let mi = gScene.materials.getMaterialInstance(sd, lod);

//    PathVertex vertex = PathVertex(0, sd.posW, sd.frame.N, sd.frontFacing);

//    int emissiveTriangleCount = gScene.lightCollection.getActiveTriangleCount();
//    if (emissiveTriangleCount == 0) return float3(0);

//    const int numEmissiveSamplesPerTriangle = 1;
//    const int numEmissiveSamples = max(1, numEmissiveSamplesPerTriangle * emissiveTriangleCount);

//    float3 direct = float3(0);

//    for (int s = 0; s < numEmissiveSamples; ++s)
//    {
//        LightSample ls;
//        if (generateEmissiveSample(vertex, true, sg, ls) && ls.pdf > 0.f)
//        {
//            SceneRayQuery<false> sceneRayQuery;
//            if (sceneRayQuery.traceVisibilityRay(ls.getVisibilityRay()))
//            {
//                direct += mi.eval(sd, ls.dir, sg)*ls.Li;
//            }
//        }
//    }

//    return direct / float(numEmissiveSamples);
//}

float3 evalDynamicDirect(ShadingData sd)
{
    let lod = ExplicitLodTextureSampler(0.f);
    let mi = gScene.materials.getMaterialInstance(sd, lod);

    PathVertex vertex = PathVertex(0, sd.posW, sd.frame.N, sd.frontFacing);

    int emissiveTriangleCount = gScene.lightCollection.getActiveTriangleCount();
    if (emissiveTriangleCount == 0) return float3(0);

    // 4 deterministic sample points per triangle.
    static const float3 kBary[4] =
    {
        float3(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0), // centroid
                                     float3(0.6, 0.2, 0.2),
        float3(0.2, 0.6, 0.2),
        float3(0.2, 0.2, 0.6)
    };

    const int kSamplesPerTri = 4;

    float3 direct = float3(0);

    for (int t = 0; t < emissiveTriangleCount; ++t)
    {
        EmissiveTriangle tri = gScene.lightCollection.getTriangle(t);

        float3 lightNormal = tri.normal;
        float triArea = tri.area;

        for (int s = 0; s < kSamplesPerTri; ++s)
        {
            float3 bary = kBary[s];

            float3 lightPos = tri.getPosition(bary);
            float2 uv = tri.getTexCoord(bary);

            float3 approxToLight = lightPos - vertex.pos;
            float3 rayOrigin = vertex.getRayOrigin(approxToLight);
            float3 toLight = lightPos - rayOrigin;

            float dist2 = dot(toLight, toLight);
            if (dist2 < 1e-8f) continue;

            float dist = sqrt(dist2);
            float3 wi = toLight / dist;

            float cosLight = max(dot(lightNormal, -wi), 0.0f);
            if (cosLight <= 0.0f) continue;

            SceneRayQuery<false> sceneRayQuery;
            if (!sceneRayQuery.traceVisibilityRay(Ray(rayOrigin, wi, 0.0f, dist - 1e-3f)))
                continue;

            float3 Le = gScene.materials.evalEmissive(tri.materialID, uv);

            // Average the fixed samples over the triangle.
            float3 Li = Le * ((triArea / float(kSamplesPerTri)) * cosLight / dist2);

            UniformSampleGenerator sg = UniformSampleGenerator(uint2(0, 0), 0);
            direct += mi.eval(sd, wi, sg) * Li;
        }
    }

    return direct;
}

struct VisibilityRayData
{
    bool visible;
};

struct PSOut
{
    float4 color : SV_Target0;
    float4 mask : SV_Target1;
};

uint hash_u32(uint x)
{
    x ^= x >> 16;
    x *= 0x7feb352d;
    x ^= x >> 15;
    x *= 0x846ca68b;
    x ^= x >> 16;
    return x;
}

PSOut psMain(VSOut vsOut, uint triangleIndex : SV_PrimitiveID)
{
    PSOut o;

    uint iid = vsOut.instanceID.index;
    if (iid < gFirstDynamicInstanceID)
    {
        discard;
    }

    float3 viewDir = normalize(gScene.camera.getPosition() - vsOut.posW);
    ShadingData shadingData = prepareShadingData(vsOut, triangleIndex, viewDir);

    BasicMaterialData baseColorMat = gScene.materials.getBasicMaterialData(vsOut.materialID);

    float4 emissiveSample = float4(baseColorMat.emissive, 1.0f);
    if (baseColorMat.texEmissive.getMode() == TextureHandle::Mode::Texture)
    {
        emissiveSample = gScene.materials.sampleTexture(baseColorMat.texEmissive, gLinearSampler, vsOut.texC, 0.f);
    }

    uint px = (uint)floor(vsOut.posW.x * 128.0f);
    uint py = (uint)floor(vsOut.posW.y * 128.0f);
    uint pz = (uint)floor(vsOut.posW.z * 128.0f);

    uint base = px * 73856093u ^ py * 19349663u ^ pz * 83492791u ^ vsOut.instanceID.index * 2654435761u;

    uint2 seed = uint2(
    hash_u32(base),
    hash_u32(base ^ 0x9e3779b9u)
);

    // UniformSampleGenerator sg = UniformSampleGenerator(uint2(vsOut.posH.xy), sampleIndex);
    UniformSampleGenerator sg = UniformSampleGenerator(seed, 0);

    float3 litColor = evalDynamicDirect(shadingData) + emissiveSample.rgb;

    o.color = float4(litColor, 1.0f);
    o.mask  = float4(1.0f, 0.0f, 0.0f, 0.0f);

    return o;
}
