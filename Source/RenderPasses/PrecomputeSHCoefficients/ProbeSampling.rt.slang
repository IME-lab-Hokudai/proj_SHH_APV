// Probe adapter for Falcor's PathTracer. The grid owns the directions and
// first-hit geometry and lighting sample budget; Falcor provides transport helpers.
import Scene.Raytracing;
import Scene.RaytracingInline;
import RenderPasses.PathTracer.PathTracer;
import RenderPasses.PathTracer.LoadShadingData;
import ProbePathTracer;
import ProbeSamplingData;
import Rendering.Lights.EmissiveLightSamplerHelpers;
import Rendering.Lights.EmissiveLightSampler;
import Utils.Geometry.GeometryHelpers;

ParameterBlock<PathTracer> gPathTracer;
StructuredBuffer<ProbeDirSample> gProbeDirSamples;
StructuredBuffer<float3> gProbePositions;
RWStructuredBuffer<ProbeSampleData> gProbeSamplingOutput;

cbuffer PerFrameCB
{
    uint probeSamplingSeed;
    uint numSamplePerProbe;      // direction samples (path traced), indices [0, numSamplePerProbe)
    uint numLightSamplePerProbe; // explicit emissive light samples, indices [numSamplePerProbe, numSamplePerProbe + numLightSamplePerProbe)
    uint numPathsPerDirection;   // K continuations per direction sharing the first hit (1 = single path, original behaviour)
    uint numNeeSamplesPerHit;    // NEE light samples at every surface hit along a path (1 = Falcor default)
};

// ---------------------------------------------------------------------------
// Clear transmissive materials (glass, bottles, liquids, ice) are invisible to probe
// rays. Probes store low-frequency diffuse irradiance: small refractive props
// contribute almost nothing to it, but traced as real glass they produce glints
// and refracted bulbs that only the fixed directions can find (spatial
// aliasing), act black after a diffuse bounce with caustics disabled, and block
// every shadow ray. Skipping them in BOTH the path rays and the visibility rays
// keeps NEE and BSDF sampling consistent (same geometry for MIS).
// Emissive materials are never skipped, even if transmissive.
// PROBE_IGNORE_TRANSMISSIVE = 0 restores the original behaviour.
// ---------------------------------------------------------------------------
#ifndef PROBE_IGNORE_TRANSMISSIVE
#define PROBE_IGNORE_TRANSMISSIVE 1
#endif

static const uint kMaxInvisibleSkips = 32;

bool isProbeInvisible(const HitInfo hit)
{
#if PROBE_IGNORE_TRANSMISSIVE
    if (hit.getType() != HitType::Triangle) return false;
    const TriangleHit th = hit.getTriangleHit();
    const uint materialID = gScene.getMaterialID(th.instanceID);
    if (gScene.materials.getMaterialType(materialID) != MaterialType::Standard) return false;
    const BasicMaterialData md = gScene.materials.getBasicMaterialData(materialID);
    const bool isEmissive = any(md.emissive > 0.f) || md.texEmissive.getMode() == TextureHandle::Mode::Texture;
    if (isEmissive) return false;
    // Clear (specular) transmission only. Diffuse-transmissive materials (thin
    // fabric, paper, leaves) scatter and block light and stay physical.
    return float(md.specularTransmission) > 0.f;
#else
    return false;
#endif
}

// Closest hit along the ray, ignoring invisible (transmissive) surfaces.
HitInfo traceSkippingInvisible(inout SceneRayQuery<kUseAlphaTest> query, Ray ray, out float hitT, uint rayFlags)
{
    hitT = 0.f;
    for (uint i = 0; i <= kMaxInvisibleSkips; ++i)
    {
        const HitInfo hit = query.traceRay(ray, hitT, rayFlags);
        if (!hit.isValid() || !isProbeInvisible(hit)) return hit;
        // Continue the same ray just past the skipped surface.
        ray.tMin = hitT * (1.f + 1e-5f) + 1e-4f;
        if (ray.tMin >= ray.tMax) break;
    }
    HitInfo miss = {};
    return miss;
}

// Visibility that lets light pass through invisible (transmissive) surfaces.
bool traceVisibilitySkippingInvisible(const Ray ray)
{
    SceneRayQuery<kUseAlphaTest> query;
    if (query.traceVisibilityRay(ray)) return true; // fast path: nothing in the way
#if PROBE_IGNORE_TRANSMISSIVE
    float hitT;
    const HitInfo hit = traceSkippingInvisible(query, ray, hitT, RAY_FLAG_NONE);
    return !hit.isValid();
#else
    return false;
#endif
}

struct ProbeVisibilityQuery : IVisibilityQuery
{
    bool traceVisibilityRay(const Ray ray)
    {
        return traceVisibilitySkippingInvisible(ray);
    }
};

// ---------------------------------------------------------------------------
// Direct light at the probe via explicit emissive light sampling.
//
// The probe position is the path origin, not a scattering vertex, so Falcor's
// NEE never runs there. Without this, direct light from small emitters is only
// found when one of the fixed directions happens to hit the emitter, which
// aliases spatially (on/off per probe). Here each probe instead takes
// numLightSamplePerProbe light samples. The sample generator is seeded by the
// sample index only (common random numbers), and the Uniform/Power emissive
// samplers pick triangles/points independently of the probe position, so every
// probe integrates over the SAME set of points on the lights. The result varies
// smoothly with probe position except at true visibility changes.
//
// Output record (consumed by accumulateDirectLightSamples on the CPU):
//   Li.rgb = emitted radiance Le (0 if occluded), Li.w = dA = 1 / (M * pdfArea)
//   s, n   = light point / light normal (polar Z-up, same .zxy convention)
//   hitT   = distance if visible, -1 if invalid or occluded (skipped on CPU)
// The CPU forms Omega_i = dA * cosXi / r^2, i.e. exactly the patch model of
// Eq. (3), so SH, gradient and Hessian use the same formulas as the ray samples.
// ---------------------------------------------------------------------------
void writeProbeLightSample(const uint lightSampleIdx, const uint outIdx, const float3 probePosW)
{
    ProbeSampleData o;
    o.Li = float4(0.f);
    o.s = float4(0.f);
    o.n = float4(0.f);
    o.hitT = -1.f;

#if USE_EMISSIVE_LIGHTS
    // Same stream for every probe: depends on the sample index only.
    SampleGenerator sg = SampleGenerator(uint2(numSamplePerProbe + lightSampleIdx, 1), probeSamplingSeed);

    TriangleLightSample tls;
    // The probe has no normal; Uniform/Power samplers ignore normal/hemisphere.
    if (gPathTracer.emissiveSampler.sampleLight(probePosW, float3(0.f, 1.f, 0.f), true, sg, tls) && tls.pdf > 0.f)
    {
        const float cosLight = dot(tls.normalW, -tls.dir);
        const float dist2 = tls.distance * tls.distance;
        if (cosLight > 1e-6f && dist2 > 1e-12f)
        {
            // Solid-angle pdf -> area pdf (position independent for Uniform/Power).
            const float pdfArea = tls.pdf * cosLight / dist2;
            const float dA = 1.f / (float(numLightSamplePerProbe) * pdfArea);

            const float3 lightPos = computeRayOrigin(tls.posW, tls.normalW);
            const float3 toLight = lightPos - probePosW;
            const float visDist = length(toLight);

            const bool visible = visDist > 0.f && traceVisibilitySkippingInvisible(Ray(probePosW, toLight / visDist, 0.f, visDist));

            o.s = float4(tls.posW.zxy, 1.f);
            o.n = float4(tls.normalW.zxy, 1.f);
            if (visible)
            {
                o.Li = float4(tls.Le, dA);
                o.hitT = tls.distance;
            }
        }
    }
#endif

    gProbeSamplingOutput[outIdx] = o;
}

[shader("raygeneration")]
void rayGen()
{
    const uint dirIdx = DispatchRaysIndex().x;
    const uint probeIdx = DispatchRaysIndex().y;
    const uint strideSamples = numSamplePerProbe + numLightSamplePerProbe;
    const uint outIdx = probeIdx * strideSamples + dirIdx;

    // Threads past the direction samples take explicit light samples.
    if (dirIdx >= numSamplePerProbe)
    {
        writeProbeLightSample(dirIdx - numSamplePerProbe, outIdx, gProbePositions[probeIdx]);
        return;
    }

    // With explicit light samples, emission seen directly by the probe (primary
    // hit) is already covered by them and must not be counted again.
    const bool includePrimaryEmission = (numLightSamplePerProbe == 0);

    // NEE light samples per surface hit. Must stay fixed along a path so the
    // emitter-hit MIS weights match the NEE at the previous vertex.
#if USE_EMISSIVE_LIGHTS
    //const uint lightSampleCount = max(1u, gScene.lightCollection.getActiveTriangleCount());
    const uint lightSampleCount = max(1u, numNeeSamplesPerHit);
#else
    const uint lightSampleCount = 1;
#endif

    // K paths per direction. All K start with the SAME primary ray, so they share
    // the first hit (s_i, n_i) that the SH gradient/Hessian patch model uses; only
    // what happens after that hit (NEE, BSDF continuation) is resampled. The
    // stored L_i is their average, i.e. a lower-variance estimate of the radiance
    // leaving patch i toward the probe. K = 1 reproduces the single-path result.
    const uint pathCount = max(1u, numPathsPerDirection);

    float3 firstHitPos = float3(0.f);
    float3 firstHitNormal = float3(0.f);
    float firstHitT = -1.f;
    float3 Lsum = float3(0.f);
    SceneRayQuery<kUseAlphaTest> sceneQuery;
    ProbeVisibilityQuery visibility;

    for (uint k = 0; k < pathCount; ++k)
    {
        PathState path = {};
        path.setActive();
        // Pack the linear DIRECTION-sample index into Falcor's 12-bit X/Y pixel fields
        // (light samples do not use PathState, so they don't count toward the limit).
        // This also supports batches with more than 4096 probes.
        path.id = probeIdx * numSamplePerProbe + dirIdx;
        path.thp = float3(1.f);
        path.origin = gProbePositions[probeIdx];
        // Preserve the grid's polar Z-up -> Falcor Y-up permutation.
        path.dir = gProbeDirSamples[dirIdx].dir.yzx;
        // Common random numbers across ALL probes, independent of batch boundaries
        // and probe ordering. Each direction (and each continuation k) initializes
        // a distinct stream; k = 0 is exactly the original stream.
        // (Light-sample threads use y = 1 with x >= numSamplePerProbe, so no overlap.)
        path.sg = SampleGenerator(uint2(dirIdx, k == 0 ? 0u : k + 1u), probeSamplingSeed);

        bool firstRay = true;

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

            // Transmissive surfaces are skipped (see isProbeInvisible), so the
            // first hit s_i is the first non-transmissive surface.
            const HitInfo hit = traceSkippingInvisible(sceneQuery, ray, hitT, rayFlags);
            if (hit.isValid())
            {
                // The primary ray is identical for all k; record its hit once.
                if (firstRay && k == 0)
                {
                    const ShadingData sd = loadShadingData(hit, path.origin, path.dir);
                    firstHitPos = sd.posW;
                    firstHitNormal = normalize(sd.faceN);
                    firstHitT = hitT;
                }
                path.setHit(hit);
                gPathTracer.handleProbeHit(path, visibility, lightSampleCount, includePrimaryEmission);
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

        Lsum += path.L;
    }

    const float3 L = Lsum / float(pathCount);

    // Preserve ProbeSampleData, output indexing, miss sentinel, and the inverse
    // coordinate permutation used by the existing SH/gradient/Hessian code.
    gProbeSamplingOutput[outIdx].Li = float4(L, 1.f);
    gProbeSamplingOutput[outIdx].s = float4(firstHitPos.zxy, 1.f);
    gProbeSamplingOutput[outIdx].n = float4(firstHitNormal.zxy, 1.f);
    gProbeSamplingOutput[outIdx].hitT = firstHitT;
}
