// FDVerification.cpp
// Finite-difference checks of the analytic SH gradient / Hessian, run from the UI.
// The maths lives in FDVerificationCore.h; this file only traces probes and writes CSV files.
//
// Four checks (FDCheck), output in the working directory:
//   FixedPatch  fd_fixedpatch.csv  random probes in the grid bounds: analytic vs fixed-patch FD (h = 1 mm)
//   Retraced    fd_fixedpatch.csv + fd_retraced.csv
//                                  also re-traced FD (h = 4 cm) and fixed-patch FD at the same step;
//                                  re-traced - fixed-patch FD isolates the visibility change
//   Plane       fd_plane_<area>.csv  every enabled area of FDPlaneSettings.h (fdplane::kRes per side):
//                                  analytic vs fixed-patch FD for all 9 coefficients, in world axes.
//                                  AdaptiveSHDemo's "Show FD areas" draws the probes in the scene.
//   PlaneRetraced fd_plane_retraced_<area>.csv  the same areas at fdplane::kRetraceRes per side:
//                                  analytic, fixed-patch FD (h = 1 mm), re-traced FD and fixed-patch FD
//                                  at kPlaneRetraceStep, plus light-visibility switches per probe.
//                                  The re-traced FD is also split by sample type: rtL re-traces only the
//                                  light samples (direction patches fixed), rtD only the direction rays
//                                  (light visibility fixed); rt - fzh = (rtL - fzh) + (rtD - fzh).
//   Both are plotted by plot_fd_plane.py.
//
// The analytic side is the pipeline itself (calculateSHCoeffsGradientsRGBAndHessiansLum and
// accumulateDirectLightSamples, single precision). The FD side is double precision, so the
// agreement floor is the single-precision rounding of the analytic side (about 1e-5).
#include "PrecomputeSHCoefficients.h"
#include "envMap_SH.h"
#include "FDVerificationCore.h"
#include "FDPlaneSettings.h"

#include <algorithm>
#include <atomic>
#include <cstdio>
#include <execution>
#include <iomanip>
#include <fstream>
#include <numeric>
#include <random>
#include <string>

namespace
{
// ---- settings --------------------------------------------------------------
// Random probes (FixedPatch, Retraced)
const uint32_t kFDProbeCount = 256;   // random probe positions inside the bounds
const uint32_t kFDPositionSeed = 1234;
const float kFDBoundsInset = 0.05f;   // keep probes 5% away from the bounds
// Fixed-patch FD verifies the implementation: no ray can switch, so a small step is safe.
const double kFixedPatchStep = 1e-3;
// Re-traced FD: rays switching objects add ~ dL (4 pi / N) / h^2 each, so the step must be large.
// The fixed-patch FD is evaluated at the same step so that the truncation error cancels.
// Scene units. The switching noise relative to the Hessian goes like (s / h)^1.5 / sqrt(N), with s the
// distance over which the lighting changes (~2 m in the SceneForData room) and N the direction count:
// with N = 4096, h of about s / 12 (here ~15 cm) gives the same noise as the 2D slice did.
const double kRetracedStep = 0.15;
// Plane checks: areas and resolution live in FDPlaneSettings.h (shared with AdaptiveSHDemo's
// "Show FD areas" view). Re-trace step of the re-traced plane check:
const double kPlaneRetraceStep = kRetracedStep;

const uint32_t kProbesPerDispatch = 256; // path-ID limit: samplesPerProbe * probes <= 2^24
// Re-traced FD: probes per chunk; a chunk holds 18 shifted traces per probe on the CPU.
const uint32_t kRetraceChunk = 16;

// Pipeline polar frame (Z up) = (world.z, world.x, world.y); world axis i is polar axis kWorldToPolar[i].
const int kWorldToPolar[3] = { 1, 2, 0 };
fdv::Vec3d toPolar(const float3& w) { return { w.z, w.x, w.y }; }
float3 toWorld(const fdv::Vec3d& p) { return float3(float(p.y), float(p.z), float(p.x)); }

void convertSamples(const ProbeSampleData* src, uint32_t count, bool lightSamples, std::vector<fdv::Sample>& dst)
{
    dst.resize(count);
    for (uint32_t i = 0; i < count; ++i)
    {
        const ProbeSampleData& s = src[i];
        fdv::Sample& d = dst[i];
        d.L = { s.Li.x, s.Li.y, s.Li.z };
        d.s = { s.s.x, s.s.y, s.s.z };
        d.n = { s.n.x, s.n.y, s.n.z };
        d.hit = s.hitT >= 0.0f;
        d.dA = lightSamples ? double(s.Li.w) : 0.0;
        if (!lightSamples && !d.hit && (s.Li.x != 0.f || s.Li.y != 0.f || s.Li.z != 0.f))
        {
            static std::atomic<bool> warned{ false };
            if (!warned.exchange(true))
                logWarning("[FD] ray miss with nonzero radiance: environment light is on; misses are treated as points at infinity.");
        }
    }
}

// Analytic derivatives exactly as the grid build computes them (float), converted to double.
// Also returns the pipeline's coefficients (luminance) for the coefficient-match check.
fdv::AnalyticDerivs analyticDerivs(const float3& xPolar, const ProbeSampleData* samples, uint32_t dirCount,
                                   uint32_t lightCount, const std::vector<ProbeDirSample>& dirs,
                                   std::array<double, 9>* coeffLum = nullptr)
{
    std::vector<float3> coeffs;
    std::vector<GradSHCoeff> grads;
    std::vector<float3x3> hess;
    calculateSHCoeffsGradientsRGBAndHessiansLum(grads, hess, xPolar, samples, dirCount, dirs);
    calculateSHCoeffs(coeffs, samples, dirCount);
    accumulateDirectLightSamples(coeffs, &grads, &hess, xPolar, samples + dirCount, lightCount);

    fdv::AnalyticDerivs an;
    for (int b = 0; b < 9; ++b)
    {
        for (int j = 0; j < 3; ++j)
        {
            an.gradRGB[b][j] = { grads[b].r[j], grads[b].g[j], grads[b].b[j] };
            for (int k = 0; k < 3; ++k) an.hessLum[b][j][k] = hess[b][j][k];
        }
        if (coeffLum) (*coeffLum)[b] = fdv::dot(fdv::Vec3d(coeffs[b].x, coeffs[b].y, coeffs[b].z), fdv::kLuma);
    }
    return an;
}

// Largest difference between the pipeline's coefficients and the FD code's coefficients at the
// probe itself, relative to the largest coefficient. Should be ~1e-6 (float): a larger value means
// the FD code evaluates a different function (basis convention, frame, sample weights).
double coefficientMatch(const std::array<double, 9>& pipelineLum, const std::vector<fdv::Sample>& dirS,
                        const std::vector<fdv::Vec3d>& dirsD, const std::vector<fdv::Sample>& lightS, const fdv::Vec3d& x0)
{
    fdv::CoeffsRGB c;
    fdv::tracedCoeffs(dirS.data(), (uint32_t)dirS.size(), dirsD.data(), lightS.data(), (uint32_t)lightS.size(), x0, c);
    double maxDiff = 0, maxAbs = 0;
    for (int b = 0; b < 9; ++b)
    {
        maxDiff = std::max(maxDiff, std::fabs(fdv::dot(c[b], fdv::kLuma) - pipelineLum[b]));
        maxAbs = std::max(maxAbs, std::fabs(pipelineLum[b]));
    }
    return maxAbs > 0 ? maxDiff / maxAbs : 0;
}

// An FD result used as the reference of compare() (for re-traced vs fixed-patch at the same step).
fdv::AnalyticDerivs asReference(const fdv::FDDerivs& d)
{
    fdv::AnalyticDerivs r;
    r.gradRGB = d.gradRGB;
    r.hessLum = d.hessLum;
    return r;
}

struct ProbeInfo
{
    double rMin = 0, hitFrac = 0;
};

ProbeInfo probeInfo(const std::vector<fdv::Sample>& dirS, const fdv::Vec3d& x0)
{
    ProbeInfo pi;
    pi.rMin = 1e30;
    uint32_t hits = 0;
    for (const auto& s : dirS)
    {
        if (!s.hit) continue;
        ++hits;
        const fdv::Vec3d q = s.s - x0;
        pi.rMin = std::min(pi.rMin, std::sqrt(fdv::dot(q, q)));
    }
    pi.hitFrac = dirS.empty() ? 0 : double(hits) / double(dirS.size());
    return pi; // rMin stays 1e30 when no direction sample hits anything
}

void writeErrorColumns(std::ofstream& f, const fdv::ErrorStats& e)
{
    f << e.gradNormAn << "," << e.gradErr << "," << e.gradRel() << ","
      << e.hessNormAn << "," << e.hessErr << "," << e.hessRel() << ","
      << e.hessRelBand(0) << "," << e.hessRelBand(1) << "," << e.hessRelBand(2);
}

// Luminance gradient (3) and Hessian (6) of coefficient b in world axes: x, y, z, xx, yy, zz, xy, xz, yz.
template<class D>
std::array<double, 9> worldComponents(const D& d, int b)
{
    std::array<double, 9> o{};
    for (int i = 0; i < 3; ++i) o[i] = fdv::dot(d.gradRGB[b][kWorldToPolar[i]], fdv::kLuma);
    const int pairs[6][2] = { { 0, 0 }, { 1, 1 }, { 2, 2 }, { 0, 1 }, { 0, 2 }, { 1, 2 } };
    for (int p = 0; p < 6; ++p) o[3 + p] = d.hessLum[b][kWorldToPolar[pairs[p][0]]][kWorldToPolar[pairs[p][1]]];
    return o;
}
} // namespace

void PrecomputeSHCoefficients::runFDVerification(RenderContext* pRenderContext, const AABB& worldBounds, FDCheck check)
{
    if (check == FDCheck::None) return;
    const uint32_t dirCount = (uint32_t)samplingDirs.size();
    const uint32_t lightCount = getLightSamplesPerProbe();
    const uint32_t stride = dirCount + lightCount;

    std::vector<fdv::Vec3d> dirsD(dirCount);
    for (uint32_t i = 0; i < dirCount; ++i)
        dirsD[i] = { samplingDirs[i].dir.x, samplingDirs[i].dir.y, samplingDirs[i].dir.z };

    // ---- buffers ----
    const uint32_t maxBatch = kProbesPerDispatch;
    mpProbePosBuffer = mpDevice->createStructuredBuffer(sizeof(float3), maxBatch, ResourceBindFlags::ShaderResource, MemoryType::DeviceLocal);
    mpProbeSamplingResultBuffer = mpDevice->createStructuredBuffer(sizeof(ProbeSampleData), size_t(stride) * maxBatch,
        ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess, MemoryType::DeviceLocal);

    // Trace a list of world positions in batches; returns stride samples per position.
    auto traceAll = [&](const float3* positions, size_t count, std::vector<ProbeSampleData>& out)
    {
        out.resize(size_t(stride) * count);
        for (size_t start = 0; start < count; start += maxBatch)
        {
            const uint32_t n = (uint32_t)std::min<size_t>(maxBatch, count - start);
            mpProbePosBuffer->setBlob(positions + start, 0, n * sizeof(float3));
            traceProbeBatch(pRenderContext, dirCount, n, lightCount);
            mpProbeSamplingResultBuffer->getBlob(out.data() + size_t(stride) * start, 0, size_t(stride) * n * sizeof(ProbeSampleData));
            mpDevice->wait();
        }
    };

    // =========================================================================
    // Plane checks: every enabled area of FDPlaneSettings.h, one CSV per area.
    // Columns: iu, iv (probe index along the in-plane axes), uAxis, vAxis (world axis index 0/1/2
    // of u and v), world position, rMin, coeffMatch, then per coefficient the luminance gradient and
    // Hessian in world axes for each estimate.
    // =========================================================================
    if (check == FDCheck::Plane || check == FDCheck::PlaneRetraced)
    {
        const bool retraced = check == FDCheck::PlaneRetraced;
        const uint32_t res = retraced ? fdplane::kRetraceRes : fdplane::kRes;
        const double h = kPlaneRetraceStep;
        const char* comp[9] = { "gx", "gy", "gz", "hxx", "hyy", "hzz", "hxy", "hxz", "hyz" };
        logInfo("[FD] {} plane check, {} x {} probes per area, fixed-patch h = {}{}; {} direction + {} light samples per probe.",
                retraced ? "re-traced" : "fixed-patch", res, res, kFixedPatchStep, retraced ? ", re-trace h = " + std::to_string(h) : std::string(),
                dirCount, lightCount);

        for (const fdplane::Area& area : fdplane::kAreas)
        {
            if (!area.enabled) continue;
            int ua, va;
            fdplane::planeAxes(area.normal, ua, va);
            const std::vector<float3> positions = fdplane::positions(area, res);
            logInfo("[FD] area {}: normal {}, centre ({}, {}, {}), half size {}, {} probes.", area.name, fdplane::normalName(area.normal),
                    area.centre.x, area.centre.y, area.centre.z, area.halfSize, positions.size());

            const std::string file = std::string(retraced ? "fd_plane_retraced_" : "fd_plane_") + area.name + ".csv";
            std::ofstream f(file);
            f << "iu,iv,uAxis,vAxis,worldX,worldY,worldZ,rMin,coeffMatch";
            if (retraced) f << ",h,dirFlipMean,lightFlipMean";
            for (int b = 0; b < 9; ++b)
                for (const char* src : retraced ? std::vector<const char*>{ "an", "fz", "rt", "fzh", "rtL", "rtD" } : std::vector<const char*>{ "an", "fz" })
                    for (const char* c : comp) f << "," << src << b << "_" << c;
            f << "\n";

            const auto offs = fdv::stencilOffsets(h);
            const uint32_t S = retraced ? uint32_t(fdv::kStencilSize) : 1u; // traces per probe
            const size_t chunk = retraced ? kRetraceChunk : maxBatch;
            std::vector<ProbeSampleData> traced;
            for (size_t p0 = 0; p0 < positions.size(); p0 += chunk)
            {
                const uint32_t pc = (uint32_t)std::min<size_t>(chunk, positions.size() - p0);
                // probe i occupies traces [i * S, i * S + S), centre first
                std::vector<float3> tracePos;
                tracePos.reserve(size_t(pc) * S);
                for (uint32_t i = 0; i < pc; ++i)
                {
                    const fdv::Vec3d x0 = toPolar(positions[p0 + i]);
                    for (uint32_t s = 0; s < S; ++s) tracePos.push_back(s == 0 ? positions[p0 + i] : toWorld(x0 + offs[s]));
                }
                traceAll(tracePos.data(), tracePos.size(), traced);

                std::vector<std::string> lines(pc);
                std::vector<uint32_t> local(pc);
                std::iota(local.begin(), local.end(), 0u);
                std::for_each(std::execution::par, local.begin(), local.end(), [&](uint32_t i)
                {
                    const size_t p = p0 + i;
                    const fdv::Vec3d x0 = toPolar(positions[p]);
                    const ProbeSampleData* s0 = traced.data() + size_t(stride) * (size_t(i) * S);
                    std::vector<fdv::Sample> dirS, lightS, ds, ls;
                    convertSamples(s0, dirCount, false, dirS);
                    convertSamples(s0 + dirCount, lightCount, true, lightS);
                    std::array<double, 9> pipeLum;
                    const fdv::AnalyticDerivs an = analyticDerivs(float3(float(x0.x), float(x0.y), float(x0.z)), s0, dirCount, lightCount, samplingDirs, &pipeLum);
                    const fdv::FDDerivs fz = fdv::frozenFD(dirS.data(), dirCount, dirsD.data(), lightS.data(), lightCount, x0, kFixedPatchStep);
                    const ProbeInfo pi = probeInfo(dirS, x0);

                    char buf[200];
                    std::string line;
                    line.reserve(retraced ? 12288 : 4096);
                    std::snprintf(buf, sizeof buf, "%u,%u,%d,%d,%.6f,%.6f,%.6f,%.6g,%.3g", uint32_t(p % res), uint32_t(p / res), ua, va,
                                  positions[p].x, positions[p].y, positions[p].z, pi.rMin, coefficientMatch(pipeLum, dirS, dirsD, lightS, x0));
                    line += buf;
                    if (!retraced)
                    {
                        for (int b = 0; b < 9; ++b)
                            for (const auto& v : { worldComponents(an, b), worldComponents(fz, b) })
                                for (double c : v) { std::snprintf(buf, sizeof buf, ",%.9g", c); line += buf; }
                        lines[i] = std::move(line);
                        return;
                    }

                    // re-traced FD from the 18 shifted traces, and fixed-patch FD at the same step
                    const fdv::FDDerivs fzh = fdv::frozenFD(dirS.data(), dirCount, dirsD.data(), lightS.data(), lightCount, x0, h);
                    // c: everything re-traced. cL: light samples re-traced, direction patches fixed.
                    // cD: direction rays re-traced, light visibility fixed (light points reprojected from x0).
                    std::array<fdv::CoeffsRGB, fdv::kStencilSize> c, cL, cD;
                    fdv::tracedCoeffs(dirS.data(), dirCount, dirsD.data(), lightS.data(), lightCount, x0, c[0]);
                    cL[0] = c[0];
                    cD[0] = c[0];
                    auto sum = [](const fdv::CoeffsRGB& a, const fdv::CoeffsRGB& b)
                    {
                        fdv::CoeffsRGB r;
                        for (int k = 0; k < 9; ++k) r[k] = a[k] + b[k];
                        return r;
                    };
                    double dMean = 0, lMean = 0;
                    for (uint32_t s = 1; s < S; ++s)
                    {
                        const ProbeSampleData* src = traced.data() + size_t(stride) * (size_t(i) * S + s);
                        convertSamples(src, dirCount, false, ds);
                        convertSamples(src + dirCount, lightCount, true, ls);
                        const fdv::Vec3d xs = x0 + offs[s];
                        fdv::tracedCoeffs(ds.data(), dirCount, dirsD.data(), ls.data(), lightCount, xs, c[s]);
                        fdv::CoeffsRGB dirFixed, lightFixed, dirRetraced, lightRetraced;
                        fdv::frozenCoeffs(dirS.data(), dirCount, dirsD.data(), nullptr, 0, x0, xs, dirFixed);
                        fdv::frozenCoeffs(dirS.data(), 0, dirsD.data(), lightS.data(), lightCount, x0, xs, lightFixed);
                        fdv::tracedCoeffs(ds.data(), dirCount, dirsD.data(), nullptr, 0, xs, dirRetraced);
                        fdv::tracedCoeffs(ds.data(), 0, dirsD.data(), ls.data(), lightCount, xs, lightRetraced);
                        cL[s] = sum(dirFixed, lightRetraced);
                        cD[s] = sum(dirRetraced, lightFixed);
                        const fdv::FlipCount fc = fdv::countFlips(dirS.data(), ds.data(), dirCount, lightS.data(), ls.data(), lightCount,
                                                                  x0, std::sqrt(fdv::dot(offs[s], offs[s])));
                        dMean += (fc.dirHits ? double(fc.dirFlips) / fc.dirHits : 0.0) / (S - 1);
                        lMean += (fc.lightValid ? double(fc.lightFlips) / fc.lightValid : 0.0) / (S - 1);
                    }
                    const fdv::FDDerivs rt = fdv::fdFromStencil(c, h);
                    const fdv::FDDerivs rtL = fdv::fdFromStencil(cL, h);
                    const fdv::FDDerivs rtD = fdv::fdFromStencil(cD, h);
                    std::snprintf(buf, sizeof buf, ",%.6g,%.6g,%.6g", h, dMean, lMean);
                    line += buf;
                    for (int b = 0; b < 9; ++b)
                        for (const auto& v : { worldComponents(an, b), worldComponents(fz, b), worldComponents(rt, b), worldComponents(fzh, b),
                                               worldComponents(rtL, b), worldComponents(rtD, b) })
                            for (double v1 : v) { std::snprintf(buf, sizeof buf, ",%.9g", v1); line += buf; }
                    lines[i] = std::move(line);
                });
                for (const auto& l : lines) f << l << "\n";
                if (!retraced || (p0 / chunk) % 10 == 0 || p0 + pc == positions.size())
                    logInfo("[FD] area {}: {} / {} probes done.", area.name, p0 + pc, positions.size());
            }
            logInfo("[FD] wrote {}", file);
        }
        return;
    }

    // =========================================================================
    // Random probes inside the bounds.
    // =========================================================================
    const bool retrace = check == FDCheck::Retraced;
    logInfo("[FD] {} random probes, {} direction + {} light samples per probe, re-trace: {}.",
            kFDProbeCount, dirCount, lightCount, retrace ? "yes" : "no");

    std::vector<float3> probesWorld(kFDProbeCount);
    {
        std::mt19937 rng(kFDPositionSeed);
        std::uniform_real_distribution<float> U(kFDBoundsInset, 1.0f - kFDBoundsInset);
        const float3 ext = worldBounds.maxPoint - worldBounds.minPoint;
        for (auto& p : probesWorld) p = worldBounds.minPoint + float3(U(rng), U(rng), U(rng)) * ext;
    }

    // ---- centre traces ----
    std::vector<ProbeSampleData> centre;
    traceAll(probesWorld.data(), probesWorld.size(), centre);

    std::vector<std::vector<fdv::Sample>> dirS(kFDProbeCount), lightS(kFDProbeCount);
    std::vector<fdv::AnalyticDerivs> analytic(kFDProbeCount);
    std::vector<ProbeInfo> info(kFDProbeCount);
    std::vector<double> coeffMatch(kFDProbeCount);
    std::vector<uint32_t> idx(kFDProbeCount);
    std::iota(idx.begin(), idx.end(), 0u);
    std::for_each(std::execution::par, idx.begin(), idx.end(), [&](uint32_t p)
    {
        const ProbeSampleData* s = centre.data() + size_t(stride) * p;
        const fdv::Vec3d x0 = toPolar(probesWorld[p]);
        convertSamples(s, dirCount, false, dirS[p]);
        convertSamples(s + dirCount, lightCount, true, lightS[p]);
        std::array<double, 9> pipeLum;
        analytic[p] = analyticDerivs(float3(float(x0.x), float(x0.y), float(x0.z)), s, dirCount, lightCount, samplingDirs, &pipeLum);
        coeffMatch[p] = coefficientMatch(pipeLum, dirS[p], dirsD, lightS[p], x0);
        info[p] = probeInfo(dirS[p], x0);
    });

    // ---- E1a: fixed-patch FD ----
    {
        struct Row { fdv::ErrorStats e; fdv::FrozenStats st; };
        std::vector<Row> rows(kFDProbeCount);
        std::for_each(std::execution::par, idx.begin(), idx.end(), [&](uint32_t p)
        {
            const fdv::Vec3d x0 = toPolar(probesWorld[p]);
            Row& r = rows[p];
            const fdv::FDDerivs fd = fdv::frozenFD(dirS[p].data(), dirCount, dirsD.data(),
                                                   lightS[p].data(), lightCount, x0, kFixedPatchStep, &r.st);
            r.e = fdv::compare(analytic[p], fd);
        });

        std::ofstream f("fd_fixedpatch.csv");
        f << std::setprecision(9);
        f << "probe,worldX,worldY,worldZ,rMin,hitFrac,coeffMatch,grazing,backFacing,h,"
             "gradNormAn,gradErr,gradRel,hessNormAn,hessErr,hessRel,hessRelL0,hessRelL1,hessRelL2\n";
        for (uint32_t p = 0; p < kFDProbeCount; ++p)
        {
            const Row& r = rows[p];
            f << p << "," << probesWorld[p].x << "," << probesWorld[p].y << "," << probesWorld[p].z << ","
              << info[p].rMin << "," << info[p].hitFrac << "," << coeffMatch[p] << ","
              << r.st.grazing << "," << r.st.backFacing << "," << kFixedPatchStep << ",";
            writeErrorColumns(f, r.e);
            f << "\n";
        }
        logInfo("[FD] wrote fd_fixedpatch.csv");
    }

    if (!retrace) return;

    // ---- E1c: re-traced FD, and fixed-patch FD at the same step ----
    // Columns: hessRel* = re-traced vs analytic (includes the truncation of the large step);
    //          hessRelVis = |re-traced - fixed-patch| / |analytic|  (visibility change only);
    //          hessRelTrunc = |fixed-patch - analytic| / |analytic| (truncation of the step only).
    std::ofstream f("fd_retraced.csv");
    f << std::setprecision(9);
    f << "probe,worldX,worldY,worldZ,rMin,hitFrac,h,"
         "gradNormAn,gradErr,gradRel,hessNormAn,hessErr,hessRel,hessRelL0,hessRelL1,hessRelL2,"
         "hessRelVis,hessRelTrunc,dirFlipMean,dirFlipMax,lightFlipMean,lightFlipMax\n";

    const double h = kRetracedStep;
    const auto offs = fdv::stencilOffsets(h);
    const uint32_t chunk = kRetraceChunk;
    for (uint32_t p0 = 0; p0 < kFDProbeCount; p0 += chunk)
    {
        const uint32_t pc = std::min(chunk, kFDProbeCount - p0);
        std::vector<float3> shifted;
        shifted.reserve(size_t(pc) * (fdv::kStencilSize - 1));
        for (uint32_t p = p0; p < p0 + pc; ++p)
        {
            const fdv::Vec3d x0 = toPolar(probesWorld[p]);
            for (int s = 1; s < fdv::kStencilSize; ++s) shifted.push_back(toWorld(x0 + offs[s]));
        }
        std::vector<ProbeSampleData> traced;
        traceAll(shifted.data(), shifted.size(), traced);

        std::vector<uint32_t> local(pc);
        std::iota(local.begin(), local.end(), 0u);
        struct Out { fdv::ErrorStats e; double vis = 0, trunc = 0, dMean = 0, dMax = 0, lMean = 0, lMax = 0; };
        std::vector<Out> outs(pc);
        std::for_each(std::execution::par, local.begin(), local.end(), [&](uint32_t i)
        {
            const uint32_t p = p0 + i;
            const fdv::Vec3d x0 = toPolar(probesWorld[p]);
            std::array<fdv::CoeffsRGB, fdv::kStencilSize> c;
            fdv::tracedCoeffs(dirS[p].data(), dirCount, dirsD.data(), lightS[p].data(), lightCount, x0, c[0]);
            std::vector<fdv::Sample> ds, ls;
            Out& o = outs[i];
            for (int s = 1; s < fdv::kStencilSize; ++s)
            {
                const ProbeSampleData* src = traced.data() + size_t(stride) * (size_t(i) * (fdv::kStencilSize - 1) + (s - 1));
                convertSamples(src, dirCount, false, ds);
                convertSamples(src + dirCount, lightCount, true, ls);
                fdv::tracedCoeffs(ds.data(), dirCount, dirsD.data(), ls.data(), lightCount, x0 + offs[s], c[s]);
                const fdv::FlipCount fc = fdv::countFlips(dirS[p].data(), ds.data(), dirCount,
                                                          lightS[p].data(), ls.data(), lightCount,
                                                          x0, std::sqrt(fdv::dot(offs[s], offs[s])));
                const double df = fc.dirHits ? double(fc.dirFlips) / fc.dirHits : 0.0;
                const double lf = fc.lightValid ? double(fc.lightFlips) / fc.lightValid : 0.0;
                o.dMean += df / (fdv::kStencilSize - 1);
                o.lMean += lf / (fdv::kStencilSize - 1);
                o.dMax = std::max(o.dMax, df);
                o.lMax = std::max(o.lMax, lf);
            }
            const fdv::FDDerivs rt = fdv::fdFromStencil(c, h);
            const fdv::FDDerivs fzh = fdv::frozenFD(dirS[p].data(), dirCount, dirsD.data(), lightS[p].data(), lightCount, x0, h);
            o.e = fdv::compare(analytic[p], rt);
            const double norm = o.e.hessNormAn > 0 ? o.e.hessNormAn : 1.0;
            o.vis = fdv::compare(asReference(fzh), rt).hessErr / norm;
            o.trunc = fdv::compare(analytic[p], fzh).hessErr / norm;
        });

        for (uint32_t i = 0; i < pc; ++i)
        {
            const uint32_t p = p0 + i;
            const Out& o = outs[i];
            f << p << "," << probesWorld[p].x << "," << probesWorld[p].y << "," << probesWorld[p].z << ","
              << info[p].rMin << "," << info[p].hitFrac << "," << h << ",";
            writeErrorColumns(f, o.e);
            f << "," << o.vis << "," << o.trunc << "," << o.dMean << "," << o.dMax << "," << o.lMean << "," << o.lMax << "\n";
        }
        logInfo("[FD] re-traced h = {}: probes {}..{} done.", h, p0, p0 + pc - 1);
    }
    logInfo("[FD] wrote fd_retraced.csv");
}
