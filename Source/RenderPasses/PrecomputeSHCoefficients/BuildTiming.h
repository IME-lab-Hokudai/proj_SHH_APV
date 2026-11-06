#pragma once
// Build-time breakdown for the placement comparison. The uniform and best-first builds fill the same counters and
// append one row per finished grid to <scene>_BuildTimings.csv, so the methods can be compared stage by stage.
//
// Wall times (ms) are measured on the calling thread around each stage and add up to the build time:
//   trace    GPU tracing of a batch, including the upload of positions and the read-back of the samples
//   shStage  the threaded per-probe step (SH, derivatives, direct light, depth) plus storing the results
//   score    scoring and splitting (best-first only)
//   other    the rest (setup, buffer allocation)
// CPU times (ms) are summed over the worker threads of parallelFor, per routine (thread time, not wall time):
//   cpuSh     SH projection of the direction samples
//   cpuDeriv  gradients, or gradients and Hessians (column "derivatives" says which)
//   cpuLight  direct-light samples (their SH and derivatives)
//   cpuDepth  harmonic-mean depth (best-first only)
// Every probe has the same samples, so the CPU time per probe compares the routines across methods, e.g. the cost
// of the Hessian = cpuDeriv per probe with Hessians - cpuDeriv per probe with gradients only.
// The timers do not change any result.
#include <atomic>
#include <chrono>
#include <cstdint>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <string>
#include <thread>
#include <algorithm>

namespace buildtiming
{
using Clock = std::chrono::steady_clock;
inline int64_t nsBetween(Clock::time_point a, Clock::time_point b)
{
    return std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
}

struct Stats
{
    uint64_t tracedProbes = 0, batches = 0;
    std::atomic<int64_t> cpuShNs{0}, cpuDerivNs{0}, cpuLightNs{0}, cpuDepthNs{0};
};

// Same thread count as probepar::parallelFor.
inline uint32_t workerThreads()
{
    return uint32_t(std::max<size_t>(1, std::min<size_t>(std::thread::hardware_concurrency(), 64)));
}

// Appends one row (header first for a new file). Wall times in ms; the CPU times come from stats.
inline void appendRow(const std::string& scene, const std::string& method, const std::string& budget,
                      const std::string& gridFile, uint64_t probes, bool hessians, uint32_t directions,
                      double buildMs, double traceMs, double shStageMs, double scoreMs, const Stats& s)
{
    const std::string path = scene + "_BuildTimings.csv";
    const bool exists = std::filesystem::exists(path);
    std::ofstream csv(path, std::ios::app);
    if (!csv) return;
    if (!exists)
        csv << "scene,method,budget,file,probes,tracedProbes,batches,directions,derivatives,threads,"
               "buildMs,traceMs,shStageMs,scoreMs,otherMs,cpuShMs,cpuDerivMs,cpuLightMs,cpuDepthMs,"
               "traceUsPerProbe,cpuShUsPerProbe,cpuDerivUsPerProbe,unixTime\n";
    const double cpuSh = s.cpuShNs.load() * 1e-6, cpuDeriv = s.cpuDerivNs.load() * 1e-6;
    const double cpuLight = s.cpuLightNs.load() * 1e-6, cpuDepth = s.cpuDepthNs.load() * 1e-6;
    const double perProbe = s.tracedProbes ? 1000.0 / double(s.tracedProbes) : 0.0; // ms -> us per traced probe
    csv << scene << "," << method << "," << budget << "," << gridFile << "," << probes << "," << s.tracedProbes << ","
        << s.batches << "," << directions << "," << (hessians ? "gradients+Hessians" : "gradients") << ","
        << workerThreads() << "," << std::fixed << std::setprecision(1) << buildMs << "," << traceMs << "," << shStageMs << ","
        << scoreMs << "," << (buildMs - traceMs - shStageMs - scoreMs) << "," << cpuSh << "," << cpuDeriv << ","
        << cpuLight << "," << cpuDepth << "," << std::setprecision(2) << traceMs * perProbe << "," << cpuSh * perProbe << ","
        << cpuDeriv * perProbe << "," << std::time(nullptr) << "\n";
}
} // namespace buildtiming
