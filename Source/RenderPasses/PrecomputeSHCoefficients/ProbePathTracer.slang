/***************************************************************************
 # Copyright (c) 2015-24, NVIDIA CORPORATION. All rights reserved.
 #
 # Redistribution and use in source and binary forms, with or without
 # modification, are permitted provided that the following conditions are met:
 #  * Redistributions of source code must retain the above copyright notice,
 #    this list of conditions and the following disclaimer.
 #  * Redistributions in binary form must reproduce the above copyright notice,
 #    this list of conditions and the following disclaimer in the documentation
 #    and/or other materials provided with the distribution.
 #  * Neither the name of NVIDIA CORPORATION nor the names of its contributors
 #    may be used to endorse or promote products derived from this software
 #    without specific prior written permission.
 #
 # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS "AS IS" AND ANY EXPRESS
 # OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 # OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN
 # NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 # INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 # BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 # DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 # OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 # NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 # EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 **************************************************************************/
import RenderPasses.PathTracer.PathTracer;
import RenderPasses.PathTracer.LoadShadingData;
import Rendering.Lights.EmissiveLightSamplerInterface;
import Rendering.Volumes.HomogeneousVolumeSampler;
import Rendering.Utils.PixelStats;

// Probe-local adaptation of Falcor's handleHit/handleMiss orchestration.
// Material evaluation, light sampling, ray offsets, dielectric state, BSDF
// continuation, and bounce limits still use the original Falcor helpers.
// Camera guide outputs and RTXDI are not part of this offline probe pass.
#if USE_RTXDI || OUTPUT_GUIDE_DATA || OUTPUT_NRD_DATA
#error ProbePathTracer requires RTXDI and camera guide outputs to be disabled.
#endif

// Diffuse-only probe transport, as in production light baking: at every surface
// hit only the diffuse lobes (reflection and transmission) are active. For the
// standard material the diffuse albedo is ~BaseColor * (1 - metallic), so metals
// come out dark; sharp specular transport (glossy reflections of lamps, caustics)
// is excluded. Clear transmissive surfaces are skipped separately
// (PROBE_IGNORE_TRANSMISSIVE). 0 = original transport.
#ifndef PROBE_DIFFUSE_ONLY
#define PROBE_DIFFUSE_ONLY 1
#endif
static const bool kProbeDiffuseOnly = PROBE_DIFFUSE_ONLY != 0;

extension PathTracer
{
    // lightSampleCount must be >= 1 and fixed along the entire path, including
    // hits rejected by nested-dielectric handling. Thus the emitter/miss MIS
    // weights use the same count as NEE at the preceding scattering vertex.
    // includePrimaryEmission: false when direct light at the probe is computed by
    // explicit light samples (ProbeSampling.rt.slang). Emission seen by the probe
    // ray itself (vertex 1) is then skipped to avoid double counting; emission
    // reached after a bounce is kept (MIS with NEE as before).
    void handleProbeHit<VisibilityQuery : IVisibilityQuery>(
        inout PathState path, inout VisibilityQuery vq, uint lightSampleCount, bool includePrimaryEmission)
    {
        // Original behaviour: Falcor's camera-style handleHit (always adds primary emission).
        if (lightSampleCount == 1 && includePrimaryEmission && !kProbeDiffuseOnly)
        {
            handleHit(path, vq);
            return;
        }

        const bool isPrimaryHit = path.getVertexIndex() == 1;
        const bool isTriangleHit = path.hit.getType() == HitType::Triangle;
        const bool isCurveHit = kUseCurves && path.hit.getType() == HitType::Curve;
        ShadingData sd = loadShadingData(path.hit, path.origin, path.dir);
        const bool isHairMaterial = kUseHairMaterial && sd.mtl.getMaterialType() == MaterialType::Hair;
        const bool isCurvePolyTubeHit = isTriangleHit && isHairMaterial;

        if (!path.interiorList.isEmpty())
        {
            const uint materialID = path.interiorList.getTopMaterialID();
            let material = gScene.materials.getMaterial(materialID);
            let vp = material.getHomogeneousVolumeProperties(gScene.materials, materialID);
            updatePathThroughput(path, HomogeneousVolumeSampler::evalTransmittance(vp, length(sd.posW - path.origin)));
        }
        if (!handleNestedDielectrics(sd, path)) return;
        logPathVertex();

        let lod = createTextureSampler(path, isPrimaryHit);
        let hints = getMaterialInstanceHints(path.hit, isPrimaryHit);
        const IMaterialInstance mi = gScene.materials.getMaterialInstance(sd, lod, hints);
        BSDFProperties bsdfProperties = mi.getProperties(sd);
        const bool isSpecular = bsdfProperties.roughness <= params.specularRoughnessThreshold;
        // Diffuse-only probe transport at every hit (including the first one),
        // or Falcor's caustic suppression after a diffuse bounce otherwise.
        if (kProbeDiffuseOnly || (kDisableCaustics && path.getBounces(BounceType::Diffuse) > 0 && isSpecular))
            sd.mtl.setActiveLobes((uint)LobeType::Diffuse);
        if (!kUseLightsInDielectricVolumes && path.isInsideDielectricVolume())
            bsdfProperties.emission = float3(0.f);

        // Emission belongs to this hit once; it is not averaged with NEE.
        const bool isLightSamplable = path.isLightSamplable();
        const bool computeEmissive = isPrimaryHit
            ? includePrimaryEmission
            : (kUseEmissiveLights && (!kUseNEE || kUseMIS || !path.isLightSampled() || !isLightSamplable));
        if (computeEmissive && any(bsdfProperties.emission > 0.f))
        {
            float misWeight = 1.f;
            if (kUseEmissiveLights && kUseNEE && kUseMIS && isTriangleHit &&
                !isPrimaryHit && path.isLightSampled() && isLightSamplable)
            {
                TriangleLightHit hit;
                TriangleHit triangleHit = path.hit.getTriangleHit();
                hit.triangleIndex = gScene.lightCollection.getTriangleIndex(triangleHit.instanceID, triangleHit.primitiveIndex);
                hit.posW = sd.posW;
                hit.normalW = sd.getOrientedFaceNormal();
                const bool upperHemisphere = path.isLightSampledUpper() && !path.isLightSampledLower();
                float lightPdf = getEmissiveSelectionProbability() *
                    emissiveSampler.evalPdf(path.origin, path.normal, upperHemisphere, hit);
                misWeight = evalMIS(1, path.pdf, float(lightSampleCount), lightPdf);
            }
            addToPathContribution(path, misWeight * bsdfProperties.emission);
        }

        // Process emitter hits at the terminal vertex before terminating,
        // since they complement the previous vertex's MIS-weighted NEE.
        if (hasFinishedSurfaceBounces(path))
        {
            path.terminate();
            return;
        }

        if (isCurveHit) path.origin = sd.posW - sd.frame.N * sd.curveRadius;
        else if (isCurvePolyTubeHit) path.origin = sd.posW + sd.frame.N * sd.curveRadius * 0.1f;
        else path.origin = sd.computeRayOrigin();

        const uint lobeTypes = mi.getLobeTypes(sd);
        const bool hasNonDeltaLobes = (lobeTypes & (uint)LobeType::NonDelta) != 0;
        path.setLightSampled(false, false);
        if (kUseNEE && hasNonDeltaLobes)
        {
            PathVertex vertex = PathVertex(path.getVertexIndex(), sd.posW, sd.faceN, sd.frontFacing);
            const bool canNeeSampleSphere = isCurveHit || isCurvePolyTubeHit;
            bool sampleUpperHemisphere = canNeeSampleSphere || ((lobeTypes & (uint)LobeType::NonDeltaReflection) != 0);
            bool sampleLowerHemisphere = canNeeSampleSphere || ((lobeTypes & (uint)LobeType::NonDeltaTransmission) != 0);
            if (!kUseLightsInDielectricVolumes && path.isInsideDielectricVolume()) sampleUpperHemisphere = false;
            path.setLightSampled(sampleUpperHemisphere, sampleLowerHemisphere);

            const float invLightSampleCount = 1.f / float(lightSampleCount);
            for (uint sampleIdx = 0; sampleIdx < lightSampleCount; ++sampleIdx)
            {
                LightSample ls = {};
                if (!generateLightSample(vertex, sampleUpperHemisphere, sampleLowerHemisphere, path.sg, ls)) continue;

                // generateLightSample already includes the light-type selection
                // probability in ls.pdf and divides ls.Li by that probability.
                if (kUseMIS && ls.lightType != (uint)LightType::Analytic)
                {
                    float scatterPdf = mi.evalPdf(sd, ls.dir, kUseBSDFSampling);
                    ls.Li *= evalMIS(float(lightSampleCount), ls.pdf, 1, scatterPdf);
                }
                float3 Lr = mi.eval(sd, ls.dir, path.sg) * ls.Li * invLightSampleCount;
                if (!any(Lr > 0.f)) continue;

                Ray ray = ls.getVisibilityRay();
                if (isCurvePolyTubeHit && dot(sd.frame.N, ray.dir) < 0.f)
                    ray.origin = ray.origin - sd.frame.N * sd.curveRadius * 2.1f;
                logTraceRay(PixelStatsRayType::Visibility);
                if (vq.traceVisibilityRay(ray)) addToPathContribution(path, Lr);
            }
            // Rejected and occluded samples still count in the denominator.
        }

        if (kUseRussianRoulette)
        {
            if (terminatePathByRussianRoulette(path, sampleNext1D(path.sg))) return;
        }

        // One BSDF sample/continuation only, after all direct-light samples.
        bool valid = generateScatterRay(sd, mi, path);
        const bool isLastVertex = hasFinishedSurfaceBounces(path);
        if (kUseNEE && !kUseMIS && isLastVertex && path.isLightSamplable()) valid = false;
        if (kDisableCaustics && path.getBounces(BounceType::Diffuse) > 0 && path.isSpecular()) valid = false;
        if (!valid) path.terminate();
    }

    void handleProbeMiss(inout PathState path, uint lightSampleCount)
    {
        if (lightSampleCount == 1)
        {
            handleMiss(path);
            return;
        }

        const bool isLightSamplable = path.isLightSamplable();
        const bool computeEnv = kUseEnvLight && (!kUseNEE || kUseMIS || !path.isLightSampled() || !isLightSamplable);
        if (computeEnv)
        {
            logPathVertex();
            float misWeight = 1.f;
            if (kUseNEE && kUseMIS && path.isLightSampled() && isLightSamplable)
            {
                float lightPdf = getEnvMapSelectionProbability() * envMapSampler.evalPdf(path.dir);
                misWeight = evalMIS(1, path.pdf, float(lightSampleCount), lightPdf);
            }
            addToPathContribution(path, misWeight * envMapSampler.eval(path.dir));
        }
        path.terminate();
    }
}
