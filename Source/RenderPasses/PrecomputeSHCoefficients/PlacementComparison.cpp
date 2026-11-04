// PlacementComparison.cpp
// Placement comparison (Comparison plan: Hessian vs. adaptive probe criteria).
//
// BestFirstBuild  every method splits its highest-scoring leaves first (AdaptiveProbeVolume::scoreCell)
//                 and stops at fixed probe budgets; a grid snapshot is saved at each budget, with one row
//                 per snapshot in <scene>_BestFirstSnapshots.csv (probe count, leaves, deepest level, the
//                 implied threshold = score of the last split, and build time split into stages).
// exportPlacementComparison
//                 reconstructs irradiance from every grid listed in <scene>_PlacementComparisonGrids.csv at a uniform
//                 lattice of points x 14 normals and compares it with the reference grid; one summary row
//                 per grid and, optionally, the per-point error. File names start with the scene name.
//
// Settings: PlacementComparisonSettings.h.
#include "PrecomputeSHCoefficients.h"
#include "PlacementComparisonSettings.h"
#include "envMap_SH.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
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

// Runs fn(i) for i in [0, n) on all hardware threads.
template<typename F>
void parallelFor(size_t n, F fn)
{
    const size_t threads = std::max<size_t>(1, std::min<size_t>(std::thread::hardware_concurrency(), 64));
    std::vector<std::thread> pool;
    const size_t chunk = (n + threads - 1) / threads;
    for (size_t t = 0; t < threads; ++t)
    {
        const size_t b = t * chunk, e = std::min(n, b + chunk);
        if (b >= e) break;
        pool.emplace_back([=, &fn]() { for (size_t i = b; i < e; ++i) fn(i); });
    }
    for (auto& th : pool) th.join();
}
} // namespace

// =====================================================================================================
// Best-first build
// =====================================================================================================
void PrecomputeSHCoefficients::BestFirstBuild(RenderContext* pRenderContext)
{
    using namespace placement;
    const auto buildStart = Clock::now();
    double traceMs = 0.0, shMs = 0.0, scoreMs = 0.0, ioMs = 0.0;
    uint32_t rounds = 0;

    std::vector<uint64_t> budgets = kProbeBudgets;
    std::sort(budgets.begin(), budgets.end());
    if (budgets.empty())
    {
        logError("[Best-first] No probe budgets set (placement::kProbeBudgets).");
        return;
    }

    auto& vol = mAdaptiveProbeVolume;
    const std::string scene = sceneTag(mpScene);
    const std::string logCsv = scene + "_BestFirstSnapshots.csv";
    const std::string scoreName = AdaptiveProbeVolume::splitScoreName(kScore);
    vol->setSplitScore(kScore);
    vol->setResidualCorrection(false);

    AABB manualBounds;
    const bool manual = getManualGridBounds(manualBounds);
    vol->startBuildSeeded(mpScene, kSeedResolution, std::numeric_limits<float>::max(), false, manual ? &manualBounds : nullptr);

    logInfo("[Best-first] Scene {}: score {}, budgets {} to {} probes, batch {}% of the leaves (at most {}).",
        scene, scoreName, budgets.front(), budgets.back(), 100.0f * kBatchFraction, kMaxBatch);

    const uint32_t dirCount = (uint32_t)samplingDirs.size();
    const uint32_t lightCount = getLightSamplesPerProbe();
    const uint32_t stride = dirCount + lightCount;
    const uint32_t kMaxCornersPerDispatch = 2048;

    mpProbePosBuffer = mpDevice->createStructuredBuffer(sizeof(float3), kMaxCornersPerDispatch, ResourceBindFlags::ShaderResource, MemoryType::DeviceLocal);
    mpProbeSamplingResultBuffer = mpDevice->createStructuredBuffer(sizeof(ProbeSampleData), size_t(stride) * kMaxCornersPerDispatch,
        ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess, MemoryType::DeviceLocal);
    std::vector<ProbeSampleData> samples(size_t(stride) * kMaxCornersPerDispatch);

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

            t0 = Clock::now();
            std::vector<std::vector<float3>> coeffs(n);
            std::vector<std::vector<GradSHCoeff>> grads(n);
            std::vector<std::vector<float3x3>> hessians(n);
            std::vector<float> hm(n);
            parallelFor(n, [&](size_t i)
            {
                const ProbeSampleData* s = samples.data() + i * size_t(stride);
                const float3 xPolar(positions[i].z, positions[i].x, positions[i].y);
                calculateSHCoeffsGradientsRGBAndHessiansLum(grads[i], hessians[i], xPolar, s, dirCount, samplingDirs);
                calculateSHCoeffs(coeffs[i], s, dirCount);
                accumulateDirectLightSamples(coeffs[i], &grads[i], &hessians[i], xPolar, s + dirCount, lightCount);
                hm[i] = harmonicMeanDepth(s, dirCount);
            });
            vol->setCornerDataRange(start, coeffs, grads, hessians);
            vol->setCornerHarmonicMeanDepthRange(start, hm);
            shMs += msSince(t0);
        }
    };

    auto saveSnapshot = [&](uint64_t budget, float impliedThreshold, bool exhausted)
    {
        const double buildMs = msSince(buildStart) - ioMs; // file I/O does not count
        const auto t0 = Clock::now();
        std::ostringstream name;
        name << scene << "_BF_" << scoreName << "_" << budget << ".txt";
        vol->setBuildTimeMs(buildMs);
        vol->saveToFile(name.str());

        const bool exists = std::ifstream(logCsv).good();
        std::ofstream csv(logCsv, std::ios::app);
        if (!exists)
            csv << "scene,score,budget,probes,leaves,deepestLevel,impliedThreshold,exhausted,rounds,"
                   "buildMs,traceMs,shMs,scoreMs,file\n";
        csv << scene << "," << scoreName << "," << budget << "," << vol->getCornerCount() << ","
            << vol->getLeafCount() << "," << vol->getDeepestLeafLevel() << "," << std::setprecision(8) << impliedThreshold << ","
            << (exhausted ? 1 : 0) << "," << rounds << "," << std::fixed << std::setprecision(1)
            << buildMs << "," << traceMs << "," << shMs << "," << scoreMs << "," << name.str() << "\n";
        ioMs += msSince(t0);
        logInfo("[Best-first] Snapshot {}: {} probes, {} leaves, deepest level {}, implied threshold {:.4g}, build {:.1f}s.",
            name.str(), vol->getCornerCount(), vol->getLeafCount(), vol->getDeepestLeafLevel(), impliedThreshold, buildMs / 1000.0);
    };

    size_t budgetIdx = 0;
    float lastScore = 0.0f;
    while (budgetIdx < budgets.size())
    {
        traceAllPending();
        ++rounds;

        const auto t0 = Clock::now();
        const auto r = vol->finishBatchBestFirst(budgets[budgetIdx], kBatchFraction, kMaxBatch);
        scoreMs += msSince(t0);
        if (r.splits > 0) lastScore = r.lastSplitScore;

        if (r.splits == 0) // budget reached or nothing left to split: the grid is complete and traced
        {
            saveSnapshot(budgets[budgetIdx], lastScore, r.exhausted);
            if (r.exhausted)
            {
                logWarning("[Best-first] Every leaf is at the maximum level; stopping before budget {}.", budgets[budgetIdx]);
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
}

// =====================================================================================================
// Evaluation
// =====================================================================================================
void PrecomputeSHCoefficients::exportPlacementComparison()
{
    using namespace placement;
    const std::string scene = sceneTag(mpScene);
    const std::string kManifestFile = scene + "_PlacementComparisonGrids.csv";
    const std::string kSummaryCsv = scene + "_PlacementComparisonSummary.csv";
    const std::string kPerPointCsv = kWritePerPointCsv ? scene + "_PlacementComparisonPoints.csv" : std::string();
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
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.empty() || line[0] == '#') continue;
            std::stringstream ss(line);
            Entry e;
            std::getline(ss, e.method, ',');
            std::getline(ss, e.type, ',');
            std::getline(ss, e.file);
            if (e.method == "method") continue; // header
            entries.push_back(e);
        }
    }
    const auto refIt = std::find_if(entries.begin(), entries.end(), [](const Entry& e) { return e.type == "reference"; });
    if (refIt == entries.end())
    {
        logError("[Placement eval] The manifest has no reference line.");
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

    std::ofstream summary(kSummaryCsv);
    summary << "method,type,file,probes,leaves,deepestLevel,buildTimeMs,samples,failedPoints,"
               "meanPctErr,p50PctErr,p95PctErr,p99PctErr,mae,rmse,meanRefLum\n";

    std::vector<std::vector<float>> perPoint; // normal-averaged percentage error per grid
    std::vector<std::string> perPointNames;

    for (const Entry& e : entries)
    {
        if (e.type == "reference") continue;
        const auto t0 = Clock::now();
        ref<AdaptiveProbeVolume> adaptive;
        ref<UniformProbeVolume> uniform;
        uint64_t probes = 0, leaves = 0; int deepest = -1; double buildMs = 0.0;
        if (e.type == "adaptive")
        {
            adaptive = AdaptiveProbeVolume::create(mpDevice);
            adaptive->loadFromFile(e.file);
            probes = adaptive->getCornerCount(); leaves = adaptive->getLeafCount();
            deepest = adaptive->getDeepestLeafLevel(); buildMs = adaptive->getBuildTimeMs();
        }
        else if (e.type == "uniform")
        {
            uniform = UniformProbeVolume::create(mpDevice);
            uniform->loadFromFile(e.file);
            const uint3 d = uniform->getProbeCountDim();
            probes = uint64_t(d.x) * d.y * d.z; leaves = uint64_t(d.x - 1) * (d.y - 1) * (d.z - 1);
            buildMs = uniform->getBuildTimeMs();
        }
        else
        {
            logWarning("[Placement eval] Unknown type '{}' for {}; skipped.", e.type, e.file);
            continue;
        }

        std::vector<float> ape(pointCount * N, std::numeric_limits<float>::quiet_NaN());
        std::vector<float> absErr(pointCount * N, std::numeric_limits<float>::quiet_NaN());
        std::vector<float> pointApe(pointCount, std::numeric_limits<float>::quiet_NaN());
        parallelFor(pointCount, [&](size_t i)
        {
            if (!refValid[i]) return;
            const float3 p = pointPos(i);
            std::vector<float3> c;
            bool ok = false;
            if (adaptive)
            {
                const int leaf = adaptive->traverseOctreeCPU(p);
                if (leaf >= 0) { adaptive->interpolateHermite_CPU(leaf, p, c); ok = c.size() >= 9; }
            }
            else ok = uniform->fetchHermiteSH_CPU(p, c) && c.size() >= 9;
            if (!ok) return;
            double sum = 0.0;
            for (size_t k = 0; k < N; ++k)
            {
                const float r = refLum[i * N + k];
                const float err = std::abs(lum709(evaluateIrradiance(normals[k], c)) - r);
                absErr[i * N + k] = err;
                ape[i * N + k] = 100.0f * err / std::max(std::abs(r), eps);
                sum += ape[i * N + k];
            }
            pointApe[i] = float(sum / double(N));
        });

        std::vector<float> valid; valid.reserve(ape.size());
        double sumAbs = 0.0, sumSq = 0.0, sumRef = 0.0; size_t failed = 0;
        for (size_t i = 0; i < pointCount; ++i)
        {
            if (!refValid[i]) continue;
            if (std::isnan(pointApe[i])) { ++failed; continue; }
            for (size_t k = 0; k < N; ++k)
            {
                valid.push_back(ape[i * N + k]);
                sumAbs += absErr[i * N + k]; sumSq += double(absErr[i * N + k]) * absErr[i * N + k];
                sumRef += refLum[i * N + k];
            }
        }
        std::sort(valid.begin(), valid.end());
        const double n = double(std::max<size_t>(1, valid.size()));
        double sumApe = 0.0; for (float v : valid) sumApe += v;

        summary << e.method << "," << e.type << "," << e.file << "," << probes << "," << leaves << "," << deepest << ","
                << std::fixed << std::setprecision(1) << buildMs << "," << valid.size() << "," << failed << ","
                << std::setprecision(6) << sumApe / n << "," << percentileSorted(valid, 50) << "," << percentileSorted(valid, 95) << ","
                << percentileSorted(valid, 99) << "," << sumAbs / n << "," << std::sqrt(sumSq / n) << "," << sumRef / n << "\n";
        summary.flush();
        logInfo("[Placement eval] {} ({} probes): mean {:.3f}%, P95 {:.3f}%, RMSE {:.4g}; {} points outside the grid; {:.1f}s.",
            e.method, probes, sumApe / n, percentileSorted(valid, 95), std::sqrt(sumSq / n), failed, msSince(t0) / 1000.0);

        if (!kPerPointCsv.empty()) { perPoint.push_back(std::move(pointApe)); perPointNames.push_back(e.method); }
    }

    if (!kPerPointCsv.empty() && !perPoint.empty())
    {
        std::ofstream pts(kPerPointCsv);
        pts << "ix,iy,iz,x,y,z,refLum";
        for (const auto& m : perPointNames) pts << "," << m;
        pts << "\n" << std::setprecision(6);
        for (size_t i = 0; i < pointCount; ++i)
        {
            if (!refValid[i]) continue;
            const float3 p = pointPos(i);
            double r = 0.0; for (size_t k = 0; k < N; ++k) r += refLum[i * N + k];
            pts << (i % R) << "," << ((i / R) % R) << "," << (i / (size_t(R) * R)) << "," << p.x << "," << p.y << "," << p.z << "," << r / double(N);
            for (const auto& v : perPoint) pts << "," << v[i];
            pts << "\n";
        }
        logInfo("[Placement eval] Wrote {}.", kPerPointCsv);
    }
    logInfo("[Placement eval] Wrote {} (absolute path: {}).", kSummaryCsv, std::filesystem::absolute(kSummaryCsv).string());
}
