// PlacementComparison.cpp
// Placement comparison (Comparison plan: Hessian vs. adaptive probe criteria).
//
// BestFirstBuild  every method splits its highest-scoring leaves first (AdaptiveProbeVolume::scoreCell)
//                 and stops at fixed budgets (probe count or runtime memory); a grid snapshot is saved at each budget, with one row
//                 per snapshot in <scene>_BestFirstSnapshots.csv (probe count, leaves, deepest level, the
//                 implied threshold = score of the last split, and build time split into stages). Every snapshot
//                 and every uniform grid of the batch also gets a row in <scene>_BuildTimings.csv (BuildTiming.h:
//                 wall time per stage and CPU time per routine, e.g. gradients vs. gradients + Hessians).
// exportPlacementComparison
//                 reconstructs irradiance from every grid listed in <scene>_PlacementComparisonGrids.csv at a uniform
//                 lattice of points x 14 normals and compares it with the reference grid; one summary row
//                 per grid and, optionally, the per-point error. File names start with the scene name.
//
// Settings: PlacementComparisonSettings.h.
#include "PrecomputeSHCoefficients.h"
#include "PlacementComparisonSettings.h"
#include "envMap_SH.h"
#include "ParallelFor.h"
#include "BuildTiming.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

// Defined in AdaptiveProbeVolume.cpp: irradiance from 9 SH coefficients for a world-space normal.
float3 evaluateIrradiance(float3 n, const std::vector<float3>& shCoeffs);

namespace
{
using Clock = std::chrono::steady_clock;
double msSince(Clock::time_point t) { return std::chrono::duration<double, std::milli>(Clock::now() - t).count(); }

// Harmonic mean of the first-hit distances of the direction samples; misses contribute 0.
float harmonicMeanDepth(const ProbeSampleData* samples, uint32_t dirCount)
{
    double invSum = 0.0;
    for (uint32_t i = 0; i < dirCount; ++i)
        if (samples[i].hitT > 0.0f) invSum += 1.0 / double(samples[i].hitT);
    return invSum > 0.0 ? float(double(dirCount) / invSum) : std::numeric_limits<float>::max();
}

// Prefix of every output file: the loaded scene's file name without extension, or the override.
std::string sceneTag(const ref<Scene>& scene)
{
    if (!placement::kSceneTagOverride.empty()) return placement::kSceneTagOverride;
    if (scene)
    {
        const auto paths = scene->getImportPaths(); // the first entry is the file that was loaded
        if (!paths.empty() && !paths.front().stem().string().empty()) return paths.front().stem().string();
    }
    return "Scene";
}

float lum709(float3 c) { return dot(c, float3(0.2126f, 0.7152f, 0.0722f)); }

// File of the reference line in a grid list (<scene>_PlacementComparisonGrids.csv); empty if none.
std::string manifestReferenceFile(const std::string& manifest)
{
    std::ifstream in(manifest);
    std::string line;
    while (std::getline(in, line))
    {
        std::vector<std::string> f;
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ','))
        {
            const auto b = cell.find_first_not_of(" \t\r\"\xEF\xBB\xBF");
            const auto e = cell.find_last_not_of(" \t\r\"");
            f.push_back(b == std::string::npos ? std::string() : cell.substr(b, e - b + 1));
        }
        if (f.size() >= 3 && f[1] == "reference") return f[2];
    }
    return std::string();
}

// Bounds from the header of a uniform grid file (UNIFORM_GRID_V1, resolution, min, max), without loading the probes.
bool readUniformGridBounds(const std::string& file, float3& bmin, float3& bmax)
{
    std::ifstream in(file);
    std::string header;
    uint32_t rx = 0, ry = 0, rz = 0;
    if (!(in >> header) || header != "UNIFORM_GRID_V1") return false;
    return bool(in >> rx >> ry >> rz >> bmin.x >> bmin.y >> bmin.z >> bmax.x >> bmax.y >> bmax.z);
}

// 6 axes + 8 diagonals (world space; the evaluators swizzle internally).
std::vector<float3> fourteenNormals()
{
    std::vector<float3> n = { float3(1, 0, 0), float3(-1, 0, 0), float3(0, 1, 0), float3(0, -1, 0), float3(0, 0, 1), float3(0, 0, -1) };
    const float d = 1.0f / std::sqrt(3.0f);
    for (int sx = -1; sx <= 1; sx += 2)
        for (int sy = -1; sy <= 1; sy += 2)
            for (int sz = -1; sz <= 1; sz += 2)
                n.push_back(float3(sx * d, sy * d, sz * d));
    return n;
}

double percentileSorted(const std::vector<float>& v, double pct)
{
    if (v.empty()) return 0.0;
    const double pos = pct * 0.01 * double(v.size() - 1);
    const size_t i0 = size_t(std::floor(pos)), i1 = std::min(i0 + 1, v.size() - 1);
    const double t = pos - double(i0);
    return double(v[i0]) * (1.0 - t) + double(v[i1]) * t;
}

using probepar::parallelFor;
} // namespace

// =====================================================================================================
// Best-first build
// =====================================================================================================
// Best-first build with the score and seed from PlacementComparisonSettings.h.
void PrecomputeSHCoefficients::BestFirstBuild(RenderContext* pRenderContext)
{
    BestFirstBuildWith(pRenderContext, placement::kScore, placement::kSeedResolution);
}

void PrecomputeSHCoefficients::BestFirstBuildWith(RenderContext* pRenderContext, AdaptiveProbeVolume::SplitScore score, uint3 seedResolution,
                                                  float gdcPercentile)
{
    using namespace placement;
    const auto buildStart = Clock::now();
    double traceMs = 0.0, shMs = 0.0, scoreMs = 0.0, ioMs = 0.0;
    uint32_t rounds = 0;
    buildtiming::Stats timing; // CPU time per routine and traced probes, for <scene>_BuildTimings.csv

    // The GDC 2005 combined rule as one best-first score (max of the three tests' ratios to their thresholds), stopped by
    // the same budgets as every other method. tau_g (gradient threshold) = gdcPercentile of the seed voxels' scores.
    const bool gdc = score == AdaptiveProbeVolume::SplitScore::GDCCombined;
    if (gdc && !(gdcPercentile >= 0.0f && gdcPercentile <= 100.0f))
    {
        logError("[Best-first] GDC combined needs a gradient percentile in [0, 100] (got {}).", gdcPercentile);
        return;
    }
    const std::vector<Budget> budgets = bestFirstBudgets(); // increasing; probe counts or runtime bytes (kBudgetUnit)
    const bool budgetInBytes = kBudgetUnit == BudgetUnit::Memory;
    if (budgets.empty())
    {
        logError("[Best-first] No budgets set (placement::kProbeBudgets / kMemoryBudgetsUniform).");
        return;
    }

    if (!mAdaptiveProbeVolume) mAdaptiveProbeVolume = AdaptiveProbeVolume::create(mpDevice); // e.g. uniform-mode builds
    auto& vol = mAdaptiveProbeVolume;
    const std::string scene = sceneTag(mpScene);
    const std::string logCsv = scene + "_BestFirstSnapshots.csv";
    // Run tag = score name, plus the seed when it is not a single root voxel (e.g. Hessian_S16), so runs with
    // different seeds write different files.
    std::string scoreName = AdaptiveProbeVolume::splitScoreName(score);
    if (seedResolution.x != 1 || seedResolution.y != 1 || seedResolution.z != 1)
    {
        scoreName += "_S" + std::to_string(seedResolution.x);
        if (seedResolution.y != seedResolution.x || seedResolution.z != seedResolution.x)
            scoreName += "x" + std::to_string(seedResolution.y) + "x" + std::to_string(seedResolution.z);
    }
    if (gdc) scoreName += "_p" + std::to_string(int(std::lround(gdcPercentile))); // e.g. GDC_S8_p90
    vol->setSplitScore(score);
    vol->setResidualCorrection(false);

    AABB manualBounds;
    const bool manual = getManualGridBounds(manualBounds);
    vol->startBuildSeeded(mpScene, seedResolution, std::numeric_limits<float>::max(), false, manual ? &manualBounds : nullptr);

    logInfo("[Best-first] Scene {}: score {}, budgets {} to {} ({}), batch {}% of the leaves (at most {}).",
        scene, scoreName, budgets.front().label, budgets.back().label,
        budgetInBytes ? "runtime memory, compact octree" : "probe count", 100.0f * kBatchFraction, kMaxBatch);

    const uint32_t dirCount = (uint32_t)samplingDirs.size();
    const uint32_t lightCount = getLightSamplesPerProbe();
    const uint32_t stride = dirCount + lightCount;
    const uint32_t kMaxCornersPerDispatch = 2048;

    mpProbePosBuffer = mpDevice->createStructuredBuffer(sizeof(float3), kMaxCornersPerDispatch, ResourceBindFlags::ShaderResource, MemoryType::DeviceLocal);
    mpProbeSamplingResultBuffer = mpDevice->createStructuredBuffer(sizeof(ProbeSampleData), size_t(stride) * kMaxCornersPerDispatch,
        ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess, MemoryType::DeviceLocal);
    std::vector<ProbeSampleData> samples(size_t(stride) * kMaxCornersPerDispatch);

    // Hessians are only computed for the scores that use them; the others store SH and gradients only (the same
    // routine as the uniform build, identical gradients), which is what Hermite interpolation needs.
    const bool needHessians = score == AdaptiveProbeVolume::SplitScore::Hessian ||
                              score == AdaptiveProbeVolume::SplitScore::HessianRelative;
    logInfo("[Best-first] Per-probe derivatives: gradients{}.", needHessians ? " and Hessians" : " only (score does not use Hessians)");
    // Per-corner data of the baselines. Interpolation always uses the analytic gradients above, the same for every
    // method; the finite-difference gradient only feeds the Gradient / GDC score.
    const bool needGeometry = score == AdaptiveProbeVolume::SplitScore::Geometry || gdc;
    const bool needFDGradient = score == AdaptiveProbeVolume::SplitScore::Gradient || gdc;
    const float3 volumeExtent = vol->getVolumeExtent();
    const float finestFraction = 1.0f / float(1u << vol->getMaxLevel());
    const float3 fdStep = volumeExtent * (kGradientFDStepCells * finestFraction); // per axis, in world units
    if (needFDGradient)
        logInfo("[Best-first] Finite-difference gradients: 6 re-traced probes per corner, step ({:.4f}, {:.4f}, {:.4f}) = {} finest cells.",
                fdStep.x, fdStep.y, fdStep.z, kGradientFDStepCells);
    uint64_t fdTraced = 0;

    // Trace every pending corner and store SH, derivatives and the harmonic-mean depth.
    auto traceAllPending = [&]()
    {
        const uint32_t total = vol->getPendingCornerCount();
        for (uint32_t start = 0; start < total; start += kMaxCornersPerDispatch)
        {
            const uint32_t count = std::min(kMaxCornersPerDispatch, total - start);
            std::vector<float3> positions;
            vol->getPendingPositionsRange(start, count, positions);
            const uint32_t n = (uint32_t)positions.size();
            if (n == 0) continue;

            auto t0 = Clock::now();
            mpProbePosBuffer->setBlob(positions.data(), 0, n * sizeof(float3));
            traceProbeBatch(pRenderContext, dirCount, n, lightCount);
            mpProbeSamplingResultBuffer->getBlob(samples.data(), 0, size_t(stride) * n * sizeof(ProbeSampleData));
            mpDevice->wait();
            traceMs += msSince(t0);
            vol->recordTraceBatch(n, dirCount);
            timing.tracedProbes += n;
            ++timing.batches;

            t0 = Clock::now();
            std::vector<std::vector<float3>> coeffs(n);
            std::vector<std::vector<GradSHCoeff>> grads(n);
            std::vector<std::vector<float3x3>> hessians(n);
            std::vector<float> hm(n);
            std::vector<std::array<float, 8>> geometry(needGeometry ? n : 0);
            parallelFor(n, [&](size_t i)
            {
                using buildtiming::nsBetween;
                const ProbeSampleData* s = samples.data() + i * size_t(stride);
                const float3 xPolar(positions[i].z, positions[i].x, positions[i].y);
                const auto c0 = Clock::now();
                if (needHessians)
                    calculateSHCoeffsGradientsRGBAndHessiansLum(grads[i], hessians[i], xPolar, s, dirCount, samplingDirs);
                else
                    calculateSHCoeffsGradients(grads[i], xPolar, s, dirCount, samplingDirs);
                const auto c1 = Clock::now();
                calculateSHCoeffs(coeffs[i], s, dirCount);
                const auto c2 = Clock::now();
                accumulateDirectLightSamples(coeffs[i], &grads[i], needHessians ? &hessians[i] : nullptr, xPolar, s + dirCount, lightCount);
                const auto c3 = Clock::now();
                hm[i] = harmonicMeanDepth(s, dirCount);
                if (needGeometry)
                {
                    // Nearest first hit per octant, as max_a |q_a| / R_a (see Corner::geometryOctant). Sampling
                    // directions are in the grid's polar frame: world = (d.y, d.z, d.x).
                    std::array<float, 8>& g = geometry[i];
                    g.fill(std::numeric_limits<float>::max());
                    for (uint32_t j = 0; j < dirCount; ++j)
                    {
                        const float t = s[j].hitT;
                        if (!(t > 0.0f)) continue;
                        const float3 d = samplingDirs[j].dir;
                        const float3 q = float3(d.y, d.z, d.x) * t;
                        const int octant = (q.x >= 0.f ? 1 : 0) | (q.y >= 0.f ? 2 : 0) | (q.z >= 0.f ? 4 : 0);
                        const float v = std::max(std::abs(q.x) / volumeExtent.x, std::max(std::abs(q.y) / volumeExtent.y, std::abs(q.z) / volumeExtent.z));
                        g[octant] = std::min(g[octant], v);
                    }
                }
                const auto c4 = Clock::now();
                timing.cpuDerivNs += nsBetween(c0, c1);
                timing.cpuShNs += nsBetween(c1, c2);
                timing.cpuLightNs += nsBetween(c2, c3);
                timing.cpuDepthNs += nsBetween(c3, c4);
            });
            vol->setCornerDataRange(start, coeffs, grads, hessians);
            vol->setCornerHarmonicMeanDepthRange(start, hm);
            if (needGeometry) vol->setCornerGeometryRange(start, geometry);
            shMs += msSince(t0);

            if (needFDGradient)
            {
                // Re-trace each corner at +-h along x, y and z (6 probes) and take central differences of the 9
                // luminance SH coefficients; the score uses the Frobenius norm. The same direction and light samples
                // are used at every position (common random numbers in ProbeSampling.rt.slang).
                const uint32_t m = 6 * n;
                std::vector<float3> shifted(m);
                for (uint32_t i = 0; i < n; ++i)
                    for (uint32_t a = 0; a < 6; ++a)
                    {
                        float3 delta(0.f);
                        const float sgn = (a & 1) ? -1.0f : 1.0f;
                        if (a < 2) delta.x = sgn * fdStep.x; else if (a < 4) delta.y = sgn * fdStep.y; else delta.z = sgn * fdStep.z;
                        shifted[6 * i + a] = positions[i] + delta;
                    }
                std::vector<std::array<float, 9>> lum(m);
                for (uint32_t s0 = 0; s0 < m; s0 += kMaxCornersPerDispatch)
                {
                    const uint32_t k = std::min(kMaxCornersPerDispatch, m - s0);
                    auto tt = Clock::now();
                    mpProbePosBuffer->setBlob(shifted.data() + s0, 0, k * sizeof(float3));
                    traceProbeBatch(pRenderContext, dirCount, k, lightCount);
                    mpProbeSamplingResultBuffer->getBlob(samples.data(), 0, size_t(stride) * k * sizeof(ProbeSampleData));
                    mpDevice->wait();
                    traceMs += msSince(tt);
                    fdTraced += k;

                    tt = Clock::now();
                    parallelFor(k, [&](size_t i)
                    {
                        const auto c0 = Clock::now();
                        const ProbeSampleData* s = samples.data() + i * size_t(stride);
                        const float3 p = shifted[s0 + i];
                        const float3 xPolar(p.z, p.x, p.y);
                        std::vector<float3> c;
                        calculateSHCoeffs(c, s, dirCount);
                        accumulateDirectLightSamples(c, nullptr, nullptr, xPolar, s + dirCount, lightCount);
                        for (int b = 0; b < 9; ++b) lum[s0 + i][b] = lum709(c[b]);
                        timing.cpuShNs += buildtiming::nsBetween(c0, Clock::now());
                    });
                    shMs += msSince(tt);
                }
                std::vector<float> gradNorm(n);
                for (uint32_t i = 0; i < n; ++i)
                {
                    double sum = 0.0;
                    const float h[3] = { fdStep.x, fdStep.y, fdStep.z };
                    for (int a = 0; a < 3; ++a)
                        for (int b = 0; b < 9; ++b)
                        {
                            const double d = (double(lum[6 * i + 2 * a][b]) - double(lum[6 * i + 2 * a + 1][b])) / (2.0 * double(h[a]));
                            sum += d * d;
                        }
                    gradNorm[i] = float(std::sqrt(sum));
                }
                vol->setCornerFDGradientRange(start, gradNorm);
            }
        }
    };

    auto saveSnapshot = [&](const Budget& budget, float impliedThreshold, bool exhausted)
    {
        const double buildMs = msSince(buildStart) - ioMs; // file I/O does not count
        const auto t0 = Clock::now();
        std::ostringstream name;
        name << scene << "_BF_" << scoreName << "_" << budget.label << ".txt";
        vol->setBuildTimeMs(buildMs);
        vol->saveToFile(name.str());

        bool exists = std::ifstream(logCsv).good();
        if (exists)
        {
            // A log from before the runtime-memory columns: keep it under another name and start a new one.
            std::string header;
            { std::ifstream in(logCsv); std::getline(in, header); }
            if (header.find("runtimeBytes") == std::string::npos)
            {
                const std::string old = scene + "_BestFirstSnapshots_old.csv";
                std::error_code ec;
                std::filesystem::rename(logCsv, old, ec);
                if (!ec) { logInfo("[Best-first] Older {} renamed to {} (new columns).", logCsv, old); exists = false; }
            }
        }
        std::ofstream csv(logCsv, std::ios::app);
        if (!exists)
            csv << "scene,score,budget,probes,leaves,deepestLevel,impliedThreshold,exhausted,rounds,"
                   "buildMs,traceMs,shMs,scoreMs,file,nodes,runtimeBytes,budgetValue\n";
        csv << scene << "," << scoreName << "," << budget.label << "," << vol->getCornerCount() << ","
            << vol->getLeafCount() << "," << vol->getDeepestLeafLevel() << "," << std::setprecision(8) << impliedThreshold << ","
            << (exhausted ? 1 : 0) << "," << rounds << "," << std::fixed << std::setprecision(1)
            << buildMs << "," << traceMs << "," << shMs << "," << scoreMs << "," << name.str() << ","
            << vol->getNodeCount() << "," << vol->getRuntimeBytes() << "," << budget.value << "\n";
        buildtiming::appendRow(scene, scoreName, budget.label, name.str(), vol->getCornerCount(), needHessians, dirCount,
                               buildMs, traceMs, shMs, scoreMs, timing);
        ioMs += msSince(t0);
        logInfo("[Best-first] Snapshot {}: {} probes, {} leaves, {:.2f} MB runtime, deepest level {}, implied threshold {:.4g}, build {:.1f}s.",
            name.str(), vol->getCornerCount(), vol->getLeafCount(), vol->getRuntimeBytes() / 1e6, vol->getDeepestLeafLevel(),
            impliedThreshold, buildMs / 1000.0);
    };

    size_t budgetIdx = 0;
    float lastScore = 0.0f;
    while (budgetIdx < budgets.size())
    {
        traceAllPending();
        ++rounds;

        if (gdc && rounds == 1)
        {
            // Gradient threshold of the combined rule (the source gives no value): this percentile of the seed voxels'
            // gradient scores (edge x max corner gradient), so the test fires on the top (100 - percentile)% of them.
            std::vector<float> g = vol->leafGradientScores();
            std::sort(g.begin(), g.end());
            const float tau = g.empty() ? 0.0f : g[std::min(g.size() - 1, size_t(std::floor(gdcPercentile / 100.0f * float(g.size() - 1) + 0.5f)))];
            vol->setGdcGradientThreshold(tau);
            logInfo("[Best-first] GDC combined: gradient threshold tau_g = P{} of {} seed voxels = {:.4g}.", gdcPercentile, g.size(), tau);
        }

        const auto t0 = Clock::now();
        const auto r = vol->finishBatchBestFirst(budgets[budgetIdx].value, kBatchFraction, kMaxBatch, budgetInBytes);
        scoreMs += msSince(t0);
        if (r.splits > 0) lastScore = r.lastSplitScore;

        if (r.splits == 0) // budget reached or nothing left to split: the grid is complete and traced
        {
            // A rule that runs out of leaves before the budget (geometry or GDC when no leaf fires any more; every leaf
            // at the maximum level) saves the grid it reached as <score>_exhausted.
            const bool early = r.exhausted;
            saveSnapshot(early ? Budget{ budgets[budgetIdx].value, "exhausted" } : budgets[budgetIdx], lastScore, r.exhausted);
            if (r.exhausted)
            {
                if (early)
                    logWarning("[Best-first] No leaf left to split (the rule no longer fires, or every leaf is at the maximum level); "
                               "stopped below budget {} and saved the grid as _exhausted.", budgets[budgetIdx].label);
                break;
            }
            ++budgetIdx;
        }
        else if (rounds % 10 == 0)
        {
            logInfo("[Best-first] Round {}: {} probes, {} leaves, last score {:.4g}; elapsed {:.1f}s.",
                rounds, vol->getCornerCount(), vol->getLeafCount(), lastScore, msSince(buildStart) / 1000.0);
        }
    }
    vol->clearPendingBatch();
    logInfo("[Best-first] Done: {} rounds, {:.1f}s (trace {:.1f}s, SH {:.1f}s, scoring {:.1f}s, file I/O {:.1f}s).",
        rounds, msSince(buildStart) / 1000.0, traceMs / 1000.0, shMs / 1000.0, scoreMs / 1000.0, ioMs / 1000.0);
    if (fdTraced > 0)
        logInfo("[Best-first] Finite differences re-traced {} probes (included in the trace and SH times).", fdTraced);
}

// =====================================================================================================
// Evaluation
// =====================================================================================================
// Entry point (UI button or kEvaluateOnLoad). Needs no scene: everything is read from grid files on the CPU.
// With a scene loaded (or kSceneTagOverride set) it evaluates <scene>_PlacementComparisonGrids.csv; with no scene it
// evaluates every *_PlacementComparisonGrids.csv in the working directory.
void PrecomputeSHCoefficients::exportPlacementComparison()
{
    static const std::string kSuffix = "_PlacementComparisonGrids.csv";
    if (mpScene || !placement::kSceneTagOverride.empty())
    {
        evaluatePlacementManifest(sceneTag(mpScene));
        return;
    }
    std::vector<std::string> tags;
    std::error_code ec;
    for (const auto& entry : std::filesystem::directory_iterator(std::filesystem::current_path(), ec))
    {
        if (!entry.is_regular_file()) continue;
        const std::string name = entry.path().filename().string();
        if (name.size() > kSuffix.size() && name.compare(name.size() - kSuffix.size(), kSuffix.size(), kSuffix) == 0)
            tags.push_back(name.substr(0, name.size() - kSuffix.size()));
    }
    if (tags.empty())
    {
        logError("[Placement eval] No scene loaded and no *{} in the working directory {}.", kSuffix,
                 std::filesystem::current_path().string());
        return;
    }
    std::sort(tags.begin(), tags.end());
    for (const std::string& t : tags)
    {
        logInfo("[Placement eval] Evaluating {}{}", t, kSuffix);
        evaluatePlacementManifest(t);
    }
}

void PrecomputeSHCoefficients::evaluatePlacementManifest(const std::string& scene)
{
    using namespace placement;
    const std::string kManifestFile = scene + "_PlacementComparisonGrids.csv";
    // One summary and per-point file per interpolation of the tested grids (kEvalInterpolations): Hermite keeps the
    // original names, trilinear adds _Trilinear. Every grid is loaded once and evaluated both ways.
    const std::vector<EvalInterpolation> modes = kEvalInterpolations;
    if (modes.empty()) { logError("[Placement eval] kEvalInterpolations is empty."); return; }
    auto modeName = [](EvalInterpolation m) { return m == EvalInterpolation::Trilinear ? std::string("trilinear") : std::string("Hermite"); };
    std::vector<std::string> summaryCsv, perPointCsv;
    for (EvalInterpolation m : modes)
    {
        const std::string suffix = m == EvalInterpolation::Trilinear ? "_Trilinear" : "";
        summaryCsv.push_back(scene + "_PlacementComparisonSummary" + suffix + ".csv");
        const bool writePoints = kWritePerPointCsv && (m == EvalInterpolation::Hermite || kWritePerPointCsvTrilinear);
        perPointCsv.push_back(writePoints ? scene + "_PlacementComparisonPoints" + suffix + ".csv" : std::string());
    }
    struct Entry { std::string method, type, file; };
    std::vector<Entry> entries;
    {
        std::ifstream in(kManifestFile);
        if (!in)
        {
            logError("[Placement eval] Manifest {} not found (absolute path: {}). One line per grid: method,type,file "
                     "with type = reference | adaptive | uniform.", kManifestFile, std::filesystem::absolute(kManifestFile).string());
            return;
        }
        std::string line;
        while (std::getline(in, line))
        {
            // Tolerate files saved by Excel: UTF-8 BOM, quoted cells, trailing ",," and padding spaces.
            auto trim = [](std::string v)
            {
                if (v.size() >= 3 && v.compare(0, 3, "\xEF\xBB\xBF") == 0) v.erase(0, 3);
                const auto notJunk = [](char ch) { return ch != ' ' && ch != '\t' && ch != '\r' && ch != '"' && ch != ','; };
                const auto b = std::find_if(v.begin(), v.end(), notJunk);
                const auto e = std::find_if(v.rbegin(), v.rend(), notJunk).base();
                return b < e ? std::string(b, e) : std::string();
            };
            line = trim(line);
            if (line.empty() || line[0] == '#') continue;
            std::stringstream ss(line);
            Entry e;
            std::getline(ss, e.method, ',');
            std::getline(ss, e.type, ',');
            std::getline(ss, e.file, ',');
            e.method = trim(e.method); e.type = trim(e.type); e.file = trim(e.file);
            if (e.method == "method") continue; // header
            if (e.type.empty() || e.file.empty())
            {
                logWarning("[Placement eval] Ignored manifest line (expected method,type,file): {}", line);
                continue;
            }
            entries.push_back(e);
        }
    }
    const auto refIt = std::find_if(entries.begin(), entries.end(), [](const Entry& e) { return e.type == "reference"; });
    if (refIt == entries.end())
    {
        logError("[Placement eval] The manifest has no reference line.");
        return;
    }

    if (!std::filesystem::exists(refIt->file))
    {
        logError("[Placement eval] Reference file {} not found.", refIt->file);
        return;
    }
    ref<UniformProbeVolume> reference = UniformProbeVolume::create(mpDevice);
    reference->loadFromFile(refIt->file);
    const float3 bmin = reference->getMinPoint(), bmax = reference->getMaxPoint();
    const uint32_t R = kEvalResolution;
    const size_t pointCount = size_t(R) * R * R;
    const std::vector<float3> normals = fourteenNormals();
    const size_t N = normals.size();

    auto pointPos = [&](size_t i)
    {
        const uint32_t x = uint32_t(i % R), y = uint32_t((i / R) % R), z = uint32_t(i / (size_t(R) * R));
        return bmin + (float3(float(x), float(y), float(z)) + 0.5f) / float(R) * (bmax - bmin);
    };

    // Reference luminance per point and normal.
    std::vector<float> refLum(pointCount * N, 0.0f);
    std::vector<uint8_t> refValid(pointCount, 0);
    parallelFor(pointCount, [&](size_t i)
    {
        std::vector<float3> c;
        if (!reference->fetchHermiteSH_CPU(pointPos(i), c) || c.size() < 9) return;
        for (size_t k = 0; k < N; ++k) refLum[i * N + k] = lum709(evaluateIrradiance(normals[k], c));
        refValid[i] = 1;
    });
    double sumAbsRef = 0.0; size_t refCount = 0;
    for (size_t i = 0; i < pointCount; ++i)
        if (refValid[i])
            for (size_t k = 0; k < N; ++k) { sumAbsRef += std::abs(refLum[i * N + k]); ++refCount; }
    const float eps = kRelEpsilonFraction * float(refCount ? sumAbsRef / double(refCount) : 1.0);
    logInfo("[Placement eval] Reference {}: {}^3 points x {} normals, epsilon {:.4g}.", refIt->file, R, N, eps);

    // Sample-point mask (classifyEvalPoints): points inside objects or near surfaces. Without the file, every point
    // counts as "away" and the ...Away columns equal the all-point ones.
    const std::string kMaskFile = scene + "_EvalPointMask.csv";
    std::vector<int8_t> pointInside(pointCount, -1);
    std::vector<float> pointMinDist(pointCount, -1.0f);
    std::vector<uint8_t> away(pointCount, 1);
    bool haveMask = false;
    {
        std::ifstream mf(kMaskFile);
        if (mf)
        {
            const float3 ext = bmax - bmin;
            const float nearDist = kNearSurfaceCells * std::min(ext.x, std::min(ext.y, ext.z)) / 64.0f; // finest cell = U64
            std::string line;
            size_t rowsRead = 0, inside = 0, nearCount = 0;
            while (std::getline(mf, line))
            {
                if (line.empty() || line[0] == '#' || line[0] == 'i') continue; // comment / header
                uint32_t ix, iy, iz, hits, backs; float minDist;
                char sep;
                std::stringstream ss(line);
                if (!(ss >> ix >> sep >> iy >> sep >> iz >> sep >> hits >> sep >> backs >> sep >> minDist)) continue;
                if (ix >= R || iy >= R || iz >= R) continue;
                const size_t i = (size_t(iz) * R + iy) * R + ix;
                const bool in = hits > 0 && 2 * backs > hits;
                const bool nr = !in && minDist >= 0.0f && minDist < nearDist;
                pointInside[i] = in ? 1 : 0;
                pointMinDist[i] = minDist;
                away[i] = (in || nr) ? 0 : 1;
                inside += in; nearCount += nr; ++rowsRead;
            }
            haveMask = rowsRead == pointCount;
            if (!haveMask)
            {
                logWarning("[Placement eval] {} has {} rows, expected {} ({}^3); ignored. Re-run the classification.", kMaskFile, rowsRead, pointCount, R);
                std::fill(away.begin(), away.end(), uint8_t(1));
                std::fill(pointInside.begin(), pointInside.end(), int8_t(-1));
                std::fill(pointMinDist.begin(), pointMinDist.end(), -1.0f);
            }
            else
                logInfo("[Placement eval] Sample-point mask {}: {} inside objects ({:.1f}%), {} within {:.3g} of a surface ({:.1f}%); "
                        "the ...Away columns use the other {:.1f}%.", kMaskFile, inside, 100.0 * inside / pointCount, nearCount, nearDist,
                        100.0 * nearCount / pointCount, 100.0 * (pointCount - inside - nearCount) / pointCount);
        }
        else
            logWarning("[Placement eval] No {} (run \"classify sample points\" with the scene loaded); ...Away columns = all points.", kMaskFile);
    }

    // Summary columns: all points, points away from geometry (...Away: kNearSurfaceCells), and points outside objects
    // (...Out: every point not inside an object; all points without the mask), the set used in the paper.
    std::vector<std::ofstream> summaries;
    for (const std::string& f : summaryCsv) summaries.emplace_back(f);
    {
        // The reference itself: zero error by definition; listed for its probe count and build time.
        const uint3 d = reference->getProbeCountDim();
        double sumRef = 0.0; size_t samplesRef = 0;
        for (size_t i = 0; i < pointCount; ++i)
            if (refValid[i])
                for (size_t k = 0; k < N; ++k) { sumRef += refLum[i * N + k]; ++samplesRef; }
        std::ostringstream row;
        row << "method,type,file,probes,leaves,deepestLevel,buildTimeMs,samples,failedPoints,"
               "meanPctErr,p50PctErr,p95PctErr,p99PctErr,mae,rmse,meanRefLum,relMSE,samplesAway,maeAway,rmseAway,relMSEAway,meanRefLumAway,"
               "memoryBytes,samplesOut,maeOut,rmseOut,relMSEOut,meanPctErrOut\n";
        row << refIt->method << ",reference," << refIt->file << "," << uint64_t(d.x) * d.y * d.z << ","
            << uint64_t(d.x - 1) * (d.y - 1) * (d.z - 1) << ",-1," << std::fixed << std::setprecision(1)
            << reference->getBuildTimeMs() << "," << samplesRef << ",0," << std::setprecision(6)
            << "0,0,0,0,0,0," << sumRef / double(std::max<size_t>(1, samplesRef)) << ",0,0,0,0,0,0,"
            << AdaptiveProbeVolume::kRuntimeProbeBytes * uint64_t(d.x) * d.y * d.z << ",0,0,0,0,0\n";
        for (auto& summary : summaries) { summary << row.str(); summary.flush(); }
    }

    // Per point and grid (normal-averaged): percentage error, absolute error, squared error, and the octree level of
    // the cell the point falls in (uniform grids: log2 of the cells per axis). Lets the error be split by region,
    // e.g. where an adaptive grid is coarser or finer than the uniform grid with the same probe count.
    // relMSE per sample = (E - E_ref)^2 / max(|E_ref|, eps)^2, the squared form of the percentage error / 100.
    struct PointColumns { std::string name; std::vector<float> pct, abs, sq, rel; std::vector<int8_t> level; };
    std::vector<std::vector<PointColumns>> perPoint(modes.size());

    for (const Entry& e : entries)
    {
        if (e.type == "reference") continue;
        if (!std::filesystem::exists(e.file))
        {
            logWarning("[Placement eval] {}: file {} not found; skipped.", e.method, e.file);
            continue;
        }
        const auto t0 = Clock::now();
        ref<AdaptiveProbeVolume> adaptive;
        ref<UniformProbeVolume> uniform;
        // memoryBytes: runtime GPU memory (uniform: probes only; adaptive: probes + compact octree + seed table).
        uint64_t probes = 0, leaves = 0, memoryBytes = 0; int deepest = -1; double buildMs = 0.0;
        if (e.type == "adaptive")
        {
            adaptive = AdaptiveProbeVolume::create(mpDevice);
            adaptive->loadFromFile(e.file);
            probes = adaptive->getCornerCount(); leaves = adaptive->getLeafCount();
            deepest = adaptive->getDeepestLeafLevel(); buildMs = adaptive->getBuildTimeMs();
            memoryBytes = adaptive->getRuntimeBytes();
        }
        else if (e.type == "uniform")
        {
            uniform = UniformProbeVolume::create(mpDevice);
            uniform->loadFromFile(e.file);
            const uint3 d = uniform->getProbeCountDim();
            probes = uint64_t(d.x) * d.y * d.z; leaves = uint64_t(d.x - 1) * (d.y - 1) * (d.z - 1);
            buildMs = uniform->getBuildTimeMs();
            memoryBytes = AdaptiveProbeVolume::kRuntimeProbeBytes * probes;
        }
        else
        {
            logWarning("[Placement eval] Unknown type '{}' for {}; skipped.", e.type, e.file);
            continue;
        }

        int uniformLevel = -1;
        if (uniform)
        {
            const uint3 d = uniform->getProbeCountDim();
            uniformLevel = int(std::lround(std::log2(double(std::max(1u, d.x - 1)))));
        }
        for (size_t mi = 0; mi < modes.size(); ++mi)
        {
        const bool trilinear = modes[mi] == EvalInterpolation::Trilinear;
        std::ofstream& summary = summaries[mi];
        std::vector<float> ape(pointCount * N, std::numeric_limits<float>::quiet_NaN());
        std::vector<float> absErr(pointCount * N, std::numeric_limits<float>::quiet_NaN());
        std::vector<float> pointApe(pointCount, std::numeric_limits<float>::quiet_NaN());
        std::vector<float> pointAbs(pointCount, std::numeric_limits<float>::quiet_NaN());
        std::vector<float> pointSq(pointCount, std::numeric_limits<float>::quiet_NaN());
        std::vector<float> pointRel(pointCount, std::numeric_limits<float>::quiet_NaN());
        std::vector<int8_t> pointLevel(pointCount, -1);
        parallelFor(pointCount, [&](size_t i)
        {
            if (!refValid[i]) return;
            const float3 p = pointPos(i);
            std::vector<float3> c;
            bool ok = false;
            if (adaptive)
            {
                const int leaf = adaptive->traverseOctreeCPU(p);
                if (leaf >= 0)
                {
                    if (trilinear) adaptive->interpolateTrilinear_CPU(leaf, p, c);
                    else adaptive->interpolateHermite_CPU(leaf, p, c);
                    ok = c.size() >= 9;
                    pointLevel[i] = int8_t(adaptive->getProbes()[leaf].level);
                }
            }
            else
            {
                ok = (trilinear ? uniform->fetchTrilinearSH_CPU(p, c) : uniform->fetchHermiteSH_CPU(p, c)) && c.size() >= 9;
                pointLevel[i] = int8_t(uniformLevel);
            }
            if (!ok) return;
            double sum = 0.0, sumA = 0.0, sumS = 0.0, sumR = 0.0;
            for (size_t k = 0; k < N; ++k)
            {
                const float r = refLum[i * N + k];
                const float err = std::abs(lum709(evaluateIrradiance(normals[k], c)) - r);
                absErr[i * N + k] = err;
                ape[i * N + k] = 100.0f * err / std::max(std::abs(r), eps);
                sum += ape[i * N + k];
                sumA += err; sumS += double(err) * err;
                const double rel = double(ape[i * N + k]) / 100.0; sumR += rel * rel;
            }
            pointApe[i] = float(sum / double(N));
            pointAbs[i] = float(sumA / double(N));
            pointSq[i] = float(sumS / double(N));
            pointRel[i] = float(sumR / double(N));
        });

        std::vector<float> valid; valid.reserve(ape.size());
        double sumAbs = 0.0, sumSq = 0.0, sumRef = 0.0; size_t failed = 0;
        double aAbs = 0.0, aSq = 0.0, aRel = 0.0, aRef = 0.0; size_t aCount = 0; // points away from geometry
        double oAbs = 0.0, oSq = 0.0, oRel = 0.0, oPct = 0.0; size_t oCount = 0;  // points outside objects
        for (size_t i = 0; i < pointCount; ++i)
        {
            if (!refValid[i]) continue;
            if (std::isnan(pointApe[i])) { ++failed; continue; }
            for (size_t k = 0; k < N; ++k)
            {
                valid.push_back(ape[i * N + k]);
                sumAbs += absErr[i * N + k]; sumSq += double(absErr[i * N + k]) * absErr[i * N + k];
                sumRef += refLum[i * N + k];
                if (pointInside[i] != 1)
                {
                    const double rel = double(ape[i * N + k]) / 100.0;
                    oAbs += absErr[i * N + k]; oSq += double(absErr[i * N + k]) * absErr[i * N + k];
                    oRel += rel * rel; oPct += ape[i * N + k]; ++oCount;
                }
                if (away[i])
                {
                    const double rel = double(ape[i * N + k]) / 100.0;
                    aAbs += absErr[i * N + k]; aSq += double(absErr[i * N + k]) * absErr[i * N + k];
                    aRel += rel * rel; aRef += refLum[i * N + k]; ++aCount;
                }
            }
        }
        const double na = double(std::max<size_t>(1, aCount));
        const double no = double(std::max<size_t>(1, oCount));
        std::sort(valid.begin(), valid.end());
        const double n = double(std::max<size_t>(1, valid.size()));
        double sumApe = 0.0, sumRel = 0.0;
        for (float v : valid) { sumApe += v; sumRel += double(v / 100.0f) * double(v / 100.0f); }

        summary << e.method << "," << e.type << "," << e.file << "," << probes << "," << leaves << "," << deepest << ","
                << std::fixed << std::setprecision(1) << buildMs << "," << valid.size() << "," << failed << ","
                << std::setprecision(6) << sumApe / n << "," << percentileSorted(valid, 50) << "," << percentileSorted(valid, 95) << ","
                << percentileSorted(valid, 99) << "," << sumAbs / n << "," << std::sqrt(sumSq / n) << "," << sumRef / n << "," << sumRel / n << ","
                << aCount << "," << aAbs / na << "," << std::sqrt(aSq / na) << "," << aRel / na << "," << aRef / na << ","
                << memoryBytes << "," << oCount << "," << oAbs / no << "," << std::sqrt(oSq / no) << "," << oRel / no << ","
                << oPct / no << "\n";
        summary.flush();
        logInfo("[Placement eval] {} ({} probes, {}): RMSE {:.4g}, MAE {:.4g}, relMSE {:.4g}; outside objects RMSE {:.4g}, "
            "MAE {:.4g}, relMSE {:.4g}; {} points outside the grid; {:.1f}s.", e.method, probes, modeName(modes[mi]),
            std::sqrt(sumSq / n), sumAbs / n, sumRel / n, std::sqrt(oSq / no), oAbs / no, oRel / no, failed, msSince(t0) / 1000.0);

        if (!perPointCsv[mi].empty())
            perPoint[mi].push_back({ e.method, std::move(pointApe), std::move(pointAbs), std::move(pointSq), std::move(pointRel), std::move(pointLevel) });
        } // interpolation modes
    }

    for (size_t mi = 0; mi < modes.size(); ++mi)
    {
        if (perPointCsv[mi].empty() || perPoint[mi].empty()) continue;
        // Columns per grid: <method>_pct, <method>_abs, <method>_sq, <method>_rel, <method>_lvl. About 400 MB for 37 grids at 64^3.
        std::ofstream pts(perPointCsv[mi]);
        pts << "ix,iy,iz,x,y,z,refLum,inside,minDist";
        for (const auto& m : perPoint[mi]) pts << "," << m.name << "_pct," << m.name << "_abs," << m.name << "_sq," << m.name << "_rel," << m.name << "_lvl";
        pts << "\n" << std::defaultfloat << std::setprecision(5);
        for (size_t i = 0; i < pointCount; ++i)
        {
            if (!refValid[i]) continue;
            const float3 p = pointPos(i);
            double r = 0.0; for (size_t k = 0; k < N; ++k) r += refLum[i * N + k];
            pts << (i % R) << "," << ((i / R) % R) << "," << (i / (size_t(R) * R)) << "," << p.x << "," << p.y << "," << p.z << "," << r / double(N) << "," << int(pointInside[i]) << "," << pointMinDist[i];
            for (const auto& m : perPoint[mi])
                pts << "," << m.pct[i] << "," << m.abs[i] << "," << m.sq[i] << "," << m.rel[i] << "," << int(m.level[i]);
            pts << "\n";
        }
        logInfo("[Placement eval] Wrote {}.", perPointCsv[mi]);
    }
    for (size_t mi = 0; mi < modes.size(); ++mi)
        logInfo("[Placement eval] Wrote {} ({} interpolation; absolute path: {}).", summaryCsv[mi], modeName(modes[mi]),
                std::filesystem::absolute(summaryCsv[mi]).string());
}

// =====================================================================================================
// Batch: build every grid of the comparison in one run, then evaluate (PlacementComparisonSettings.h)
// =====================================================================================================
// Uniform grid of cellsPerAxis^3 cells, same tracing and threaded SH/gradient step as the uniform build in
// execute(); saved to file with its build time (file I/O excluded).
void PrecomputeSHCoefficients::buildUniformGridToFile(RenderContext* pRenderContext, uint32_t cellsPerAxis, const std::string& file)
{
    const auto tStart = Clock::now();
    ref<UniformProbeVolume> vol = UniformProbeVolume::create(mpDevice);
    AABB bounds;
    const bool manual = getManualGridBounds(bounds);
    vol->initGrid(mpScene, uint3(cellsPerAxis, cellsPerAxis, cellsPerAxis), manual ? &bounds : nullptr);
    const uint3 dim = vol->getProbeCountDim();
    const uint32_t totalProbes = dim.x * dim.y * dim.z;

    const uint32_t dirCount = (uint32_t)samplingDirs.size();
    const uint32_t lightCount = getLightSamplesPerProbe();
    const uint32_t stride = dirCount + lightCount;
    const uint32_t kMaxProbesPerDispatch = 2048;
    mpProbePosBuffer = mpDevice->createStructuredBuffer(sizeof(float3), kMaxProbesPerDispatch, ResourceBindFlags::ShaderResource, MemoryType::DeviceLocal);
    mpProbeSamplingResultBuffer = mpDevice->createStructuredBuffer(sizeof(ProbeSampleData), size_t(stride) * kMaxProbesPerDispatch,
        ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess, MemoryType::DeviceLocal);

    std::vector<ProbeSampleData> samples(size_t(stride) * kMaxProbesPerDispatch);
    std::vector<float3> positions(kMaxProbesPerDispatch);
    std::vector<std::vector<float3>> coeffs(kMaxProbesPerDispatch);
    std::vector<std::vector<GradSHCoeff>> grads(kMaxProbesPerDispatch);
    double traceMs = 0.0, shMs = 0.0;
    buildtiming::Stats timing; // same counters as the best-first build, for <scene>_BuildTimings.csv
    uint32_t nextProgress = 25;
    for (uint32_t batchStart = 0; batchStart < totalProbes; batchStart += kMaxProbesPerDispatch)
    {
        const uint32_t n = std::min(kMaxProbesPerDispatch, totalProbes - batchStart);
        for (uint32_t i = 0; i < n; ++i)
        {
            const uint32_t probeIdx = batchStart + i; // (z * Ny + y) * Nx + x, as in the uniform build
            const uint32_t x = probeIdx % dim.x, y = (probeIdx / dim.x) % dim.y, z = probeIdx / (dim.x * dim.y);
            positions[i] = vol->getMinPoint() + float3(x, y, z) * vol->getCellSize();
        }
        auto t0 = Clock::now();
        mpProbePosBuffer->setBlob(positions.data(), 0, n * sizeof(float3));
        traceProbeBatch(pRenderContext, dirCount, n, lightCount);
        pRenderContext->submit(true);
        mpProbeSamplingResultBuffer->getBlob(samples.data(), 0, size_t(stride) * n * sizeof(ProbeSampleData));
        mpDevice->wait();
        traceMs += msSince(t0);
        timing.tracedProbes += n;
        ++timing.batches;

        t0 = Clock::now();
        parallelFor(n, [&](size_t i)
        {
            using buildtiming::nsBetween;
            const float3 p = positions[i];
            const float3 xPolar(p.z, p.x, p.y);
            const ProbeSampleData* s = samples.data() + i * size_t(stride);
            const auto c0 = Clock::now();
            calculateSHCoeffsGradients(grads[i], xPolar, s, dirCount, samplingDirs);
            const auto c1 = Clock::now();
            calculateSHCoeffs(coeffs[i], s, dirCount);
            const auto c2 = Clock::now();
            accumulateDirectLightSamples(coeffs[i], &grads[i], nullptr, xPolar, s + dirCount, lightCount);
            const auto c3 = Clock::now();
            timing.cpuDerivNs += nsBetween(c0, c1);
            timing.cpuShNs += nsBetween(c1, c2);
            timing.cpuLightNs += nsBetween(c2, c3);
        });
        for (uint32_t i = 0; i < n; ++i) vol->setProbeData(batchStart + i, coeffs[i], grads[i]);
        shMs += msSince(t0);

        const uint32_t done = batchStart + n;
        if (100ull * done >= uint64_t(nextProgress) * totalProbes)
        {
            logInfo("[Batch] U{}: {}/{} probes, {:.1f}s.", cellsPerAxis, done, totalProbes, msSince(tStart) / 1000.0);
            nextProgress += 25;
        }
    }
    vol->setBuildTimeMs(msSince(tStart));
    vol->saveToFile(file);
    buildtiming::appendRow(sceneTag(mpScene), "Uniform", "U" + std::to_string(cellsPerAxis), file, totalProbes, false, dirCount,
                           vol->getBuildTimeMs(), traceMs, shMs, 0.0, timing);
    logInfo("[Batch] U{} done: {} probes in {:.1f}s (trace {:.1f}s, SH {:.1f}s) -> {}", cellsPerAxis, totalProbes,
            vol->getBuildTimeMs() / 1000.0, traceMs / 1000.0, shMs / 1000.0, file);
}

void PrecomputeSHCoefficients::runPlacementBatch(RenderContext* pRenderContext)
{
    using namespace placement;
    const auto t0 = Clock::now();
    const std::string scene = sceneTag(mpScene);
    std::vector<std::string> failed;
    {
        AABB b;
        if (getManualGridBounds(b))
            logInfo("[Batch] Scene {}: probe volume ({:.2f}, {:.2f}, {:.2f}) to ({:.2f}, {:.2f}, {:.2f}).", scene,
                    b.minPoint.x, b.minPoint.y, b.minPoint.z, b.maxPoint.x, b.maxPoint.y, b.maxPoint.z);
        else logInfo("[Batch] Scene {}: probe volume = scene bounds with the 2% inset.", scene);
    }
    auto step = [&](const std::string& name, const std::function<void()>& fn)
    {
        const auto ts = Clock::now();
        logInfo("[Batch] ==== {} (elapsed {:.1f} min) ====", name, msSince(t0) / 60000.0);
        try { fn(); }
        catch (const std::exception& e) { logError("[Batch] {} failed: {}", name, e.what()); failed.push_back(name); }
        catch (...) { logError("[Batch] {} failed.", name); failed.push_back(name); }
        logInfo("[Batch] ==== {} finished in {:.1f} min ====", name, msSince(ts) / 60000.0);
    };

    for (uint32_t res : kBatchUniformResolutions)
    {
        const std::string file = scene + "_U" + std::to_string(res) + ".txt";
        if (kBatchSkipExisting && std::filesystem::exists(file))
        {
            logInfo("[Batch] uniform U{}: {} exists, skipped.", res, file);
            continue;
        }
        step("uniform U" + std::to_string(res), [&]() { buildUniformGridToFile(pRenderContext, res, file); });
    }
    // Reference grid of this scene (kSceneSetups), built only when its file is missing.
    if (const SceneSetup* setup = findSceneSetup(scene); setup && setup->referenceResolution > 0)
    {
        const uint32_t res = setup->referenceResolution;
        const std::string file = scene + "_U" + std::to_string(res) + ".txt";
        if (std::filesystem::exists(file)) logInfo("[Batch] reference U{}: {} exists, skipped.", res, file);
        else step("reference U" + std::to_string(res), [&]() { buildUniformGridToFile(pRenderContext, res, file); });
    }
    for (AdaptiveProbeVolume::SplitScore sc : kBatchScores)
    {
        // Same run tag as BestFirstBuildWith: score name, plus _S<seed> for a seed other than 1 x 1 x 1.
        std::string tag = AdaptiveProbeVolume::splitScoreName(sc);
        if (kSeedResolution.x != 1 || kSeedResolution.y != 1 || kSeedResolution.z != 1)
        {
            tag += "_S" + std::to_string(kSeedResolution.x);
            if (kSeedResolution.y != kSeedResolution.x || kSeedResolution.z != kSeedResolution.x)
                tag += "x" + std::to_string(kSeedResolution.y) + "x" + std::to_string(kSeedResolution.z);
        }
        const std::vector<Budget> budgets = bestFirstBudgets();
        const std::string lastFile = scene + "_BF_" + tag + "_" + (budgets.empty() ? std::string("0") : budgets.back().label) + ".txt";
        const std::string exhaustedFile = scene + "_BF_" + tag + "_exhausted.txt"; // a rule that stopped below the last budget
        if (kBatchSkipExisting && (std::filesystem::exists(lastFile) || std::filesystem::exists(exhaustedFile)))
        {
            logInfo("[Batch] best-first {}: {} exists, skipped.", tag, std::filesystem::exists(lastFile) ? lastFile : exhaustedFile);
            continue;
        }
        step(std::string("best-first ") + AdaptiveProbeVolume::splitScoreName(sc),
             [&]() { BestFirstBuildWith(pRenderContext, sc, kSeedResolution); });
    }
    // GDC 2005 combined rule at the same budgets, one run per gradient-threshold percentile.
    for (float pct : kBatchGdcPercentiles)
    {
        std::string tag = std::string(AdaptiveProbeVolume::splitScoreName(AdaptiveProbeVolume::SplitScore::GDCCombined));
        if (kSeedResolution.x != 1 || kSeedResolution.y != 1 || kSeedResolution.z != 1)
        {
            tag += "_S" + std::to_string(kSeedResolution.x);
            if (kSeedResolution.y != kSeedResolution.x || kSeedResolution.z != kSeedResolution.x)
                tag += "x" + std::to_string(kSeedResolution.y) + "x" + std::to_string(kSeedResolution.z);
        }
        tag += "_p" + std::to_string(int(std::lround(pct)));
        const std::vector<Budget> budgets = bestFirstBudgets();
        const std::string lastFile = scene + "_BF_" + tag + "_" + (budgets.empty() ? std::string("0") : budgets.back().label) + ".txt";
        const std::string exhaustedFile = scene + "_BF_" + tag + "_exhausted.txt";
        if (kBatchSkipExisting && (std::filesystem::exists(lastFile) || std::filesystem::exists(exhaustedFile)))
        {
            logInfo("[Batch] GDC combined {}: {} exists, skipped.", tag, std::filesystem::exists(lastFile) ? lastFile : exhaustedFile);
            continue;
        }
        step("GDC combined " + tag,
             [&]() { BestFirstBuildWith(pRenderContext, AdaptiveProbeVolume::SplitScore::GDCCombined, kSeedResolution, pct); });
    }
    if (kBatchClassifyPoints)
        step("classify sample points", [&]() { classifyEvalPoints(pRenderContext); });
    if (kBatchEvaluateAtEnd)
        step("evaluation", [&]() { exportPlacementComparison(); });

    if (failed.empty()) logInfo("[Batch] All steps done in {:.1f} min.", msSince(t0) / 60000.0);
    else
    {
        std::string list; for (const auto& f : failed) list += (list.empty() ? "" : ", ") + f;
        logError("[Batch] Done in {:.1f} min; failed steps: {}", msSince(t0) / 60000.0, list);
    }
}

// =====================================================================================================
// Sample-point classification: which evaluation points are inside objects or near surfaces
// =====================================================================================================
// Same points as the evaluation (kEvalResolution^3 cell centres of the reference grid's bounds). For each point,
// kClassifyDirections rays without back-face culling (ProbeSampling.rt.slang, classifyMode = 1) give the number of
// hits, the number of back-face hits and the closest hit distance. Written to <scene>_EvalPointMask.csv; the
// evaluation turns them into "inside" / "near" flags (PlacementComparisonSettings.h).
void PrecomputeSHCoefficients::classifyEvalPoints(RenderContext* pRenderContext)
{
    using namespace placement;
    const auto t0 = Clock::now();
    if (!mpScene) { logError("[Placement classify] Needs a loaded scene."); return; }
    const std::string scene = sceneTag(mpScene);
    const std::string manifest = scene + "_PlacementComparisonGrids.csv";
    const std::string refFile = manifestReferenceFile(manifest);
    float3 bmin, bmax;
    if (refFile.empty() || !readUniformGridBounds(refFile, bmin, bmax))
    {
        logError("[Placement classify] Could not read the reference grid bounds (grid list {}, reference '{}').", manifest, refFile);
        return;
    }
    const uint32_t R = kEvalResolution;
    const size_t pointCount = size_t(R) * R * R;
    const uint32_t K = kClassifyDirections;
    const uint32_t kBatch = 2048;

    mpProbePosBuffer = mpDevice->createStructuredBuffer(sizeof(float3), kBatch, ResourceBindFlags::ShaderResource, MemoryType::DeviceLocal);
    mpProbeSamplingResultBuffer = mpDevice->createStructuredBuffer(sizeof(ProbeSampleData), size_t(K) * kBatch,
        ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess, MemoryType::DeviceLocal);
    std::vector<ProbeSampleData> samples(size_t(K) * kBatch);
    std::vector<float3> positions(kBatch);
    std::vector<uint16_t> hits(pointCount, 0), backs(pointCount, 0);
    std::vector<float> minDist(pointCount, -1.0f);

    struct ModeGuard { uint32_t& m; ~ModeGuard() { m = 0; } } guard{ mProbeClassifyMode };
    mProbeClassifyMode = 1;
    for (size_t start = 0; start < pointCount; start += kBatch)
    {
        const uint32_t n = uint32_t(std::min<size_t>(kBatch, pointCount - start));
        for (uint32_t j = 0; j < n; ++j)
        {
            const size_t i = start + j;
            const uint32_t x = uint32_t(i % R), y = uint32_t((i / R) % R), z = uint32_t(i / (size_t(R) * R));
            positions[j] = bmin + (float3(float(x), float(y), float(z)) + 0.5f) / float(R) * (bmax - bmin); // as in the evaluation
        }
        mpProbePosBuffer->setBlob(positions.data(), 0, n * sizeof(float3));
        traceProbeBatch(pRenderContext, K, n, 0);
        pRenderContext->submit(true);
        mpProbeSamplingResultBuffer->getBlob(samples.data(), 0, size_t(K) * n * sizeof(ProbeSampleData));
        mpDevice->wait();
        for (uint32_t j = 0; j < n; ++j)
        {
            const size_t i = start + j;
            for (uint32_t k = 0; k < K; ++k)
            {
                const ProbeSampleData& sd = samples[size_t(j) * K + k];
                if (sd.hitT < 0.0f) continue;
                ++hits[i];
                if (sd.Li.x > 0.5f) ++backs[i];
                if (minDist[i] < 0.0f || sd.hitT < minDist[i]) minDist[i] = sd.hitT;
            }
        }
    }
    mProbeClassifyMode = 0;

    const std::string out = scene + "_EvalPointMask.csv";
    std::ofstream f(out);
    f << "# Evaluation-point classification: " << R << "^3 points in [" << bmin.x << "," << bmin.y << "," << bmin.z << "] - ["
      << bmax.x << "," << bmax.y << "," << bmax.z << "], " << K << " rays each, no back-face culling.\n";
    f << "ix,iy,iz,hits,backHits,minDist\n" << std::setprecision(6);
    size_t inside = 0;
    for (size_t i = 0; i < pointCount; ++i)
    {
        f << (i % R) << "," << ((i / R) % R) << "," << (i / (size_t(R) * R)) << "," << hits[i] << "," << backs[i] << "," << minDist[i] << "\n";
        inside += (hits[i] > 0 && 2 * backs[i] > hits[i]);
    }
    logInfo("[Placement classify] {} points, {} rays each: {} inside objects ({:.1f}%). Wrote {} in {:.1f}s.",
        pointCount, K, inside, 100.0 * inside / pointCount, out, msSince(t0) / 1000.0);
}
