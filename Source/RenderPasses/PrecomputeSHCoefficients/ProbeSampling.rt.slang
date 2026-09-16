// Probe adapter for Falcor's PathTracer. The grid owns the directions and
// first-hit geometry and lighting sample budget; Falcor provides transport helpers.
import Scene.Raytracing;
import Scene.RaytracingInline;
import RenderPasses.PathTracer.PathTracer;
import RenderPasses.PathTracer.LoadShadingData;
import ProbePathTracer;
import ProbeSamplingData;

ParameterBlock<PathTracer> gPathTracer;
StructuredBuffer<ProbeDirSample> gProbeDirSamples;
StructuredBuffer<float3> gProbePositions;
RWStructuredBuffer<ProbeSampleData> gProbeSamplingOutput;

cbuffer PerFrameCB
{
    uint probeSamplingSeed;
    uint numSamplePerProbe;
};

struct ProbeVisibilityQuery : IVisibilityQuery
{
    bool traceVisibilityRay(const Ray ray)
    {
        SceneRayQuery<kUseAlphaTest> query;
        return query.traceVisibilityRay(ray);
    }
};

// Fixed offline budget per surface hit, independent of scene triangle count.
// 1 uses Falcor's original hit/miss handlers; larger values use local multi-sample NEE.
static const uint kProbeDirectLightSamples = 1;

[shader("raygeneration")]
void rayGen()
{
    const uint dirIdx = DispatchRaysIndex().x;
    const uint probeIdx = DispatchRaysIndex().y;
    const uint outIdx = probeIdx * numSamplePerProbe + dirIdx;

    PathState path = {};
    path.setActive();
    // Pack the linear sample index into Falcor's 12-bit X/Y pixel fields.
    // This also supports batches with more than 4096 probes.
    path.id = outIdx;
    path.thp = float3(1.f);
    path.origin = gProbePositions[probeIdx];
    // Preserve the grid's polar Z-up -> Falcor Y-up permutation.
    path.dir = gProbeDirSamples[dirIdx].dir.yzx;
    // Common random numbers across ALL probes, independent of batch boundaries
    // and probe ordering. Each direction still initializes a distinct stream.
    path.sg = SampleGenerator(uint2(dirIdx, 0), probeSamplingSeed);

    float3 firstHitPos = float3(0.f);
    float3 firstHitNormal = float3(0.f);
    float firstHitT = -1.f;
    bool firstRay = true;
    SceneRayQuery<kUseAlphaTest> sceneQuery;
    ProbeVisibilityQuery visibility;
    const uint lightSampleCount = max(1u, kProbeDirectLightSamples);

    // The probe is a point measuring radiance, not a diffuse receiver.
    // Start at vertex zero with unit throughput and no cosine/PDF weighting.
    while (path.isActive())
    {
        path.incrementVertexIndex();
        Ray ray = path.getScatterRay();
        if (firstRay)
        {
            ray.tMin = 0.001f;
            ray.tMax = 10000.f;
        }
        float hitT;
        // Preserve first-ray culling; secondary rays must be able to exit glass.
        const uint rayFlags = firstRay
            ? RAY_FLAG_CULL_BACK_FACING_TRIANGLES
            : RAY_FLAG_NONE;

        const HitInfo hit = sceneQuery.traceRay(ray, hitT, rayFlags);
        if (hit.isValid())
        {
            if (firstRay)
            {
                const ShadingData sd = loadShadingData(hit, path.origin, path.dir);
                firstHitPos = sd.posW;
                firstHitNormal = normalize(sd.faceN);
                firstHitT = hitT;
            }
            path.setHit(hit);
            gPathTracer.handleProbeHit(path, visibility, lightSampleCount);
        }
        else
        {
            path.clearHit();
            gPathTracer.handleProbeMiss(path, lightSampleCount);
        }
        // Capture only the first intersection, even when later hits are rejected
        // by Falcor's nested dielectric logic or the ray transmits through glass.
        firstRay = false;
    }

    // Preserve ProbeSampleData, output indexing, miss sentinel, and the inverse
    // coordinate permutation used by the existing SH/gradient/Hessian code.
    gProbeSamplingOutput[outIdx].Li = float4(path.L, 1.f);
    gProbeSamplingOutput[outIdx].s = float4(firstHitPos.zxy, 1.f);
    gProbeSamplingOutput[outIdx].n = float4(firstHitNormal.zxy, 1.f);
    gProbeSamplingOutput[outIdx].hitT = firstHitT;
}
