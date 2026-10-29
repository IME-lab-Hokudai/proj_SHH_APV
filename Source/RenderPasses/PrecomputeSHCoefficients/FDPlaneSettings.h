// FDPlaneSettings.h
// The areas where the FD plane checks put their probes. Shared by FDVerification.cpp
// (PrecomputeSHCoefficients) and the "Show FD areas" view in AdaptiveSHDemo, so the dots drawn
// there are exactly the probes the checks evaluate. Edit the areas here only.
//
// Each area is a square of res x res probes in an axis-aligned plane. The checks run every enabled
// area and write one CSV per area (fd_plane_<name>.csv, fd_plane_retraced_<name>.csv).
#pragma once
#include "Falcor.h"
#include <array>
#include <cstdint>
#include <vector>

namespace fdplane
{
using Falcor::float3;

// Probes per side: fixed-patch check (FDCheck::Plane) and re-traced check (FDCheck::PlaneRetraced).
// The re-traced check costs ~ res^2 * 19 probe traces per area.
inline constexpr uint32_t kRes = 100;
inline constexpr uint32_t kRetraceRes = 50;

enum class Normal { X, Y, Z }; // plane normal; Y = horizontal square

struct Area
{
    const char* name;  // used in file names
    bool enabled;
    Normal normal;
    float3 centre;     // world, metres
    float halfSize;    // half the side length
    float3 color;      // AdaptiveSHDemo dot colour
};

// Defaults for SceneForData.pyscene (room x, z in [-5, 5], y in [0, 4]); three lighting regimes.
//   A  y = 1, beside ThinSlabE (12.5 cm from it), half under the green light: penumbrae from the
//      warm, green and wall lights cross the square.
//   B  y = 3, 1 m under the warm ceiling light, nothing in between: strong, smooth curvature
//      (closest geometry: the ceiling 1 m above; TallBoxA's top is 0.4 m below and outside the square).
//   C  y = 1, behind ThinSlabF near the right wall (11 cm from the slab, 15 cm from the wall): deep
//      shadow of the cool and green lights with its penumbra edge.
// (The PG Fig. 4 region, centre (0, 0.2, 0), half size 0.25, cuts 2.5 cm into ThinSlabE in this room.)
inline const std::array<Area, 3> kAreas = { {
    { "A_slabE_penumbra", true, Normal::Y, float3(0.9f, 1.0f, -0.8f), 1.0f, float3(1.0f, 0.45f, 0.0f) },
    { "B_under_warm_light", true, Normal::Y, float3(-2.8f, 3.0f, -1.8f), 1.0f, float3(0.0f, 0.9f, 1.0f) },
    { "C_slabF_shadow", true, Normal::Y, float3(4.35f, 1.0f, 0.8f), 0.5f, float3(1.0f, 0.2f, 0.9f) },
} };

// In-plane axes (u, v) as world axis indices: Y -> (x, z), X -> (z, y), Z -> (x, y).
inline void planeAxes(Normal n, int& u, int& v)
{
    if (n == Normal::Y) { u = 0; v = 2; }
    else if (n == Normal::X) { u = 2; v = 1; }
    else { u = 0; v = 1; }
}
inline char axisName(int a) { return a == 0 ? 'x' : (a == 1 ? 'y' : 'z'); }
inline char normalName(Normal n) { return n == Normal::X ? 'x' : (n == Normal::Y ? 'y' : 'z'); }

// Probe positions (world, Y up), row by row: index = iv * res + iu.
inline std::vector<float3> positions(const Area& a, uint32_t res)
{
    int u, v;
    planeAxes(a.normal, u, v);
    const double step = 2.0 * a.halfSize / double(res - 1);
    std::vector<float3> p;
    p.reserve(size_t(res) * res);
    for (uint32_t iv = 0; iv < res; ++iv)
        for (uint32_t iu = 0; iu < res; ++iu)
        {
            float3 q = a.centre;
            q[u] = float(a.centre[u] - a.halfSize + step * iu);
            q[v] = float(a.centre[v] - a.halfSize + step * iv);
            p.push_back(q);
        }
    return p;
}
inline float spacing(const Area& a, uint32_t res) { return 2.0f * a.halfSize / float(res - 1); }
} // namespace fdplane
