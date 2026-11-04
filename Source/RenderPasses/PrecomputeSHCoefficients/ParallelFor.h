#pragma once
// Runs fn(i) for i in [0, n) on all hardware threads (at most 64), each thread taking one contiguous
// range. fn must only read shared data and write to its own index i; results are then identical to
// a serial loop. Used for the per-probe SH / derivative computation after each traced batch.
#include <algorithm>
#include <cstddef>
#include <thread>
#include <vector>

namespace probepar
{
template<typename F>
void parallelFor(size_t n, F&& fn)
{
    const size_t threads = std::max<size_t>(1, std::min<size_t>(std::thread::hardware_concurrency(), 64));
    if (n == 0) return;
    if (threads == 1 || n == 1)
    {
        for (size_t i = 0; i < n; ++i) fn(i);
        return;
    }
    const size_t chunk = (n + threads - 1) / threads;
    std::vector<std::thread> pool;
    pool.reserve(threads);
    for (size_t t = 0; t < threads; ++t)
    {
        const size_t b = t * chunk, e = std::min(n, b + chunk);
        if (b >= e) break;
        pool.emplace_back([b, e, &fn]() { for (size_t i = b; i < e; ++i) fn(i); });
    }
    for (auto& th : pool) th.join();
}
} // namespace probepar
