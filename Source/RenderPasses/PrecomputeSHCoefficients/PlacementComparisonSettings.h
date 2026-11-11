#pragma once
// Settings of the placement comparison (Comparison plan: Hessian vs. adaptive probe criteria).
// Shared by PrecomputeSHCoefficients.cpp (build switch) and PlacementComparison.cpp.
#include "AdaptiveProbeVolume.h"
#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>

namespace placement
{
// ---- Build -------------------------------------------------------------------------------------
// Threshold: the original build (finishBatch, ErrorThreshold, residual correction as configured).
// BestFirst: every method splits its highest-scoring leaves first and stops at fixed probe budgets;
//            a snapshot grid is saved at each budget (residual correction off).
enum class BuildMode { Threshold, BestFirst };
inline constexpr BuildMode kBuildMode = BuildMode::BestFirst;

// Score of the best-first build: Hessian (ours, E_abs), HessianRelative (ours, E_rel; floor = 1% of the seed's
// mean ||L||) or HarmonicMeanDepth.
//inline constexpr AdaptiveProbeVolume::SplitScore kScore = AdaptiveProbeVolume::SplitScore::HarmonicMeanDepth;
//inline constexpr AdaptiveProbeVolume::SplitScore kScore = AdaptiveProbeVolume::SplitScore::Hessian;
inline constexpr AdaptiveProbeVolume::SplitScore kScore = AdaptiveProbeVolume::SplitScore::HessianRelative;

// Seed grid of the best-first build. Levels are absolute (AdaptiveProbeVolume::mMaxLevel), so the seed does not
// change the finest cell size; 1 x 1 x 1 matches the PG paper's single root voxel. 8 x 8 x 8 starts from U8
// (729 probes) and spends the rest of each budget on refinement. A seed other than 1 x 1 x 1 is added to the
// file names: <scene>_BF_<score>_S8_<budget>.txt.
//inline const uint3 kSeedResolution = uint3(1, 1, 1);
//inline const uint3 kSeedResolution = uint3(16, 16, 16);
inline const uint3 kSeedResolution = uint3(8, 8, 8);

// What a budget limits.
//   Memory: the runtime GPU memory of the grid (AdaptiveProbeVolume::runtimeBytes: probes, octree nodes in the
//           compact 36-byte layout, seed table). Budgets are given as uniform grid resolutions (kMemoryBudgetsUniform):
//           budget n = the memory of uniform U<n> = 576 B x (n + 1)^3 probes, so every snapshot has a uniform grid of
//           equal memory. Files: <scene>_BF_<score>_S8_memU<n>.txt.
//   Probes: the probe count (kProbeBudgets). An adaptive grid also stores its octree, so equal probes is not equal
//           memory (about 5% more with the compact nodes). Files: <scene>_BF_<score>_S8_<budget>.txt.
enum class BudgetUnit { Probes, Memory };
inline constexpr BudgetUnit kBudgetUnit = BudgetUnit::Memory;

// Uniform resolutions whose memory is a budget. The U8 seed alone (729 probes + 512 nodes) is above U8's memory,
// so the smallest is U16. Each needs a uniform grid of the same resolution for the comparison (kBatchUniformResolutions).
inline const std::vector<uint32_t> kMemoryBudgetsUniform = { 16, 20, 26, 32, 40, 51 };

// Probe budgets = total probes of the grid, seed included (probe = grid corner).
// U8 = 9^3, U16 = 17^3, U32 = 33^3, U64 = 65^3. With the U8 seed, the 729 snapshot is the seed itself
// (a check: it must match a uniform U8 grid).
inline const std::vector<uint64_t> kProbeBudgets = { 729, 4913, 10000, 20000, 35937, 70000, 140000 };

// Budgets of the best-first build in increasing order, with the label used in the snapshot file names.
struct Budget { uint64_t value; std::string label; };
inline std::vector<Budget> bestFirstBudgets()
{
    std::vector<Budget> b;
    if (kBudgetUnit == BudgetUnit::Memory)
        for (uint32_t n : kMemoryBudgetsUniform)
            b.push_back({ AdaptiveProbeVolume::kRuntimeProbeBytes * uint64_t(n + 1) * (n + 1) * (n + 1), "memU" + std::to_string(n) });
    else
        for (uint64_t n : kProbeBudgets) b.push_back({ n, std::to_string(n) });
    std::sort(b.begin(), b.end(), [](const Budget& x, const Budget& y) { return x.value < y.value; });
    return b;
}

// Leaves split per round: this fraction of the current leaf count, at most kMaxBatch.
inline constexpr float kBatchFraction = 0.01f;
inline constexpr uint32_t kMaxBatch = 256;

// Every output file starts with the scene tag: the loaded scene's file name without extension
// (e.g. SceneForData, Bistro), or kSceneTagOverride when it is not empty.
//   <scene>_BF_<score>[_S<seed>]_<budget>.txt   snapshot grid at each budget
//   <scene>_BestFirstSnapshots.csv    one row per snapshot
inline const std::string kSceneTagOverride = "";

// ---- Per-scene setup ---------------------------------------------------------------------------
// Probe-volume bounds and reference grid per scene (matched by the scene tag above). Every grid of the comparison
// (uniform, best-first, reference) and the evaluation lattice use these bounds, so all cover the same region.
// A scene not listed here keeps the old behaviour: the scene bounds with the 2% inset (or kGridMin/kGridMax in
// PrecomputeSHCoefficients.cpp when kUseManualGridBounds is set there).
//   boundsMin/boundsMax  region in Falcor world coordinates (Y up), used when manualBounds is true
//   referenceResolution  > 0: the batch also builds <scene>_U<n>.txt as the reference, once (skipped when the file
//                        exists, even with kBatchSkipExisting = false, because it is the most expensive grid)
struct SceneSetup
{
    std::string scene;
    bool manualBounds;
    float3 boundsMin, boundsMax;
    uint32_t referenceResolution;
};
inline const std::vector<SceneSetup> kSceneSetups = {
    // Bistro interior (scene turned to the axes in the pyscene, rotationYDeg about (9.52, 0, -1.08)): a 6.5 x 5.35 x 6.5 m
    // part of the dining area, floor to ceiling (pendant lamps and both ceiling fans inside), centred at x 9.52,
    // z 0.42: the sphere-loop demo cube's centre moved 1.5 m in +z, off the bar counter (-z) onto the tables. Near-cubic, so
    // the cells are near-cubic too (U64: about 10 x 8.4 x 10 cm, close to SceneForData's 7.7 cm). The larger TC region
    // (3.5, 0.65, -6.5)-(16, 5, 3) gave cells 2.3x longer in x than in y, coarser than the lighting detail.
    { "BistroInterior_Wine", true, float3(6.27f, 0.65f, -2.83f), float3(12.77f, 6.0f, 3.67f), 128 },
};
inline const SceneSetup* findSceneSetup(const std::string& scene)
{
    for (const auto& s : kSceneSetups)
        if (s.scene == scene) return &s;
    return nullptr;
}

// ---- Evaluation --------------------------------------------------------------------------------
// Per scene, in Mogwai's working directory:
//   <scene>_PlacementComparisonGrids.csv     input: one grid per line, method,type,file
//                                            (type = reference | adaptive | uniform; exactly one reference,
//                                            a uniform grid such as U64 or U128; lines starting with # skipped)
//   <scene>_PlacementComparisonSummary.csv   output: one row per grid
//   <scene>_PlacementComparisonPoints.csv    output: error per point and grid (if kWritePerPointCsv)
inline constexpr bool kWritePerPointCsv = true;
// Keep the results of grids evaluated before: a grid is not evaluated again (not even loaded) when its row is in the
// existing summary of each interpolation (same columns), neither its file nor the reference is newer than that summary,
// and its columns are in the existing per-point file; the kept rows and columns are copied into the new files.
// Set false (or delete the summaries) after changing the evaluation itself, e.g. the interpolation code.
inline constexpr bool kEvalSkipEvaluated = true;

// How the tested grids are interpolated at the evaluation points (the reference is read at its own probes, where both
// give the probe value). One summary per entry:
//   Hermite:   tricubic Hermite with the analytic gradients (the runtime method); <scene>_PlacementComparisonSummary.csv
//   Trilinear: corner values only, no gradients;  <scene>_PlacementComparisonSummary_Trilinear.csv
//   LimitedHermite: Hermite with monotone (Fritsch-Carlson) slope limiting, so a cell cannot overshoot the values at
//              its corners (e.g. next to an emitter, where a probe's slope only holds over its distance to the light);
//              <scene>_PlacementComparisonSummary_LimitedHermite.csv
// Trilinear separates placement from the gradients, which leave out the visibility boundary term.
enum class EvalInterpolation { Hermite, Trilinear, LimitedHermite };
inline const std::vector<EvalInterpolation> kEvalInterpolations = { EvalInterpolation::Hermite, EvalInterpolation::Trilinear,
                                                                    EvalInterpolation::LimitedHermite };
// Octree seams. Where a leaf meets a finer leaf, plain interpolation of the coarse leaf ignores the finer probes on the
// shared face, so the lighting jumps across the face (Oat, "Irradiance Volumes for Games", GDC 2005, slides 42-43).
//   None: plain interpolation of the leaf containing the point (as before; files without _Oat).
//   Oat:  Oat's principle (slide 43: keep using the samples of finer neighbours), made exact: in a lattice cell of the
//         finest level that touches a finer leaf, the point is interpolated on that lattice cell, whose 8 nodes take
//         their value (and, for Hermite, gradient) from the smallest leaf touching each node. Neighbouring lattice cells
//         share their nodes, so faces between leaves are continuous (exactly for trilinear; for Hermite up to a small
//         residual where a lattice cell meets the leaf's own interpolation). Lattice cells whose nodes all belong to the
//         leaf keep the leaf's own interpolation. Files get _Oat; uniform grids are unchanged. Applied to every
//         interpolation in kEvalInterpolations.
// The summary reports the seam: the mean and largest jump of irradiance luminance (14 normals) across faces between
// leaves of the same or different levels (16 points per face; seam... columns), and which seam handling was used.
//   Constrain: hanging-node constraint (AdaptiveProbeVolume::constrainHangingCorners): every corner on a face or edge of a
//         coarser leaf takes that leaf's value there (and its blended gradient), so the grid itself is seamless and the
//         runtime needs nothing extra (no lookups, no memory); uses the coarse side's data on those faces. Files get
//         _Constrained; uniform grids are unchanged.
enum class SeamFix { None, Oat, Constrain };
// Oat is dropped from the comparison (2026-11-11): its runtime lookups cost the demo about 3x the frame time, too slow to
// report; the hanging-node constraint is free at runtime. None = the raw grid (no seam handling), for comparison.
// Add SeamFix::Oat here to evaluate it again.
inline const std::vector<SeamFix> kEvalSeamFixes = { SeamFix::None, SeamFix::Constrain };

// Per-point files of the other interpolations (<scene>_PlacementComparisonPoints_<mode>.csv, each as large as the
// Hermite one, about 400 MB for SceneForData).
inline constexpr bool kWritePerPointCsvTrilinear = true;

inline constexpr uint32_t kEvalResolution = 64;   // uniform lattice of cell centres per axis
inline constexpr float kRelEpsilonFraction = 0.01f; // epsilon of the percentage error, x mean reference luminance

// Run the evaluation once when the scene loads (also available as a UI button).
inline constexpr bool kEvaluateOnLoad = false;

// ---- Batch: everything in one run --------------------------------------------------------------
// Builds the uniform grids, then one best-first run per score (seed kSeedResolution, budgets kProbeBudgets),
// then evaluates <scene>_PlacementComparisonGrids.csv. Start it with the UI button
// "Placement comparison: build everything + evaluate", or set kRunBatchOnLoad to start it when the scene loads.
// Uniform grids are saved as <scene>_U<cells>.txt (e.g. SceneForData_U32.txt). Progress is in the log ([Batch]).
inline constexpr bool kRunBatchOnLoad = false;
// One uniform grid per memory budget (kMemoryBudgetsUniform), plus U8 and U64. U20, U26, U40 and U51 also sit near
// the old probe budgets 10,000, 20,000, 70,000 and 140,000.
inline const std::vector<uint32_t> kBatchUniformResolutions = { 8, 16, 20, 26, 32, 40, 51, 64 };
inline const std::vector<AdaptiveProbeVolume::SplitScore> kBatchScores = {
    AdaptiveProbeVolume::SplitScore::HarmonicMeanDepth,
    AdaptiveProbeVolume::SplitScore::Hessian,
    AdaptiveProbeVolume::SplitScore::HessianRelative,
    // Baselines (GDC 2005 / Greger et al. 1998). Geometry stops once no voxel with geometry can be split further and
    // saves that grid as <scene>_BF_Geometry_S8_exhausted.txt (plus every budget it reached before).
    AdaptiveProbeVolume::SplitScore::Geometry,
    // Oat's gradient test as written (GDC 2005, slides 21-22 and 35: "large gradient", gradients by central differences
    // of 6 extra samples): max corner finite-difference gradient norm, no size factor.
    AdaptiveProbeVolume::SplitScore::GradientRaw,
    // Our extension of it, as an additional experiment: e x max corner FD gradient, the first-order change across the
    // voxel (a size-aware version of the same test).
    AdaptiveProbeVolume::SplitScore::Gradient,
};
// GDC 2005 combined rule at the same memory budgets: best-first on max(e / HM_min, e ||grad|| / tau_g, geometry ? 1 : 0),
// so every test fires at 1 and the budget replaces that threshold. tau_g sets how strongly the gradient test counts
// against depth and geometry; the source gives no value, so one run per percentile of the seed voxels' e ||grad||
// (p50: the gradient test fires on half the seed voxels; p90: on the top 10%).
// Saved as <scene>_BF_GDC_S8_p<percentile>_memU<n>.txt. Empty = not built. Each run traces 7x the probes.
//inline const std::vector<float> kBatchGdcPercentiles = { 50, 90 };
inline const std::vector<float> kBatchGdcPercentiles = {}; // off for now: Geometry and Gradient first
// Finite-difference step of the Gradient and GDC scores, in finest-level (U64) cells along each axis. The probes are
// re-traced at +-h, so rays that switch objects add noise ~ 1/h: a step much below a cell gets noisy (FDVerification).
inline constexpr float kGradientFDStepCells = 1.0f;
inline constexpr bool kBatchEvaluateAtEnd = true;
// Skip grids whose files already exist (uniform: the grid file; best-first: the snapshot at the largest budget),
// so adding a resolution or a score only builds what is missing. Set false to rebuild everything (e.g. new timings);
// the batch then also starts new <scene>_BuildTimings.csv and <scene>_BestFirstSnapshots.csv (the old ones are kept
// as *_previous.csv), so old and new runs are not mixed.
// Rebuild of 2026-11-11 (Bistro lamp textures on the build PC, corner keys on the finest lattice): false, both scenes.
inline constexpr bool kBatchSkipExisting = false;
// The reference grid (kSceneSetups) is the most expensive one and is normally built only when its file is missing,
// whatever kBatchSkipExisting says. True: rebuild it as well (Bistro's U128 had the old lamp textures).
inline constexpr bool kBatchRebuildReference = true;
inline constexpr bool kBatchClassifyPoints = true; // classify the evaluation points before the evaluation (scene loaded)

// ---- Sample points away from geometry ------------------------------------------------------------
// classifyEvalPoints traces kClassifyDirections rays (no back-face culling) from every evaluation point and writes
// <scene>_EvalPointMask.csv. A point is "inside" when most rays that hit something hit a back face, and "near" when
// the closest hit is nearer than kNearSurfaceCells finest cells (U64 cell). The evaluation then also reports every
// metric over the points that are neither (columns ...Away in the summary; inside and minDist in the per-point file).
inline constexpr uint32_t kClassifyDirections = 64;
// One finest cell: closer than that, neighbouring probes can lie on both sides of a surface and interpolation is not
// meaningful (with 0.5 a few such points remained and skewed relMSE in coarse regions).
inline constexpr float kNearSurfaceCells = 1.0f;

// View evaluation (offline images, PlacementComparison.cpp): the evaluation also reconstructs the irradiance at the
// pixels of every <scene>_View_*.dump written by AdaptiveSHDemo ("Dump view"), for every grid and mode.
inline constexpr bool kEvalViews = true;
// Each pixel is evaluated (and its reference traced) this far off the surface along its normal, in scene units
// (metres in Bistro), so a dynamic object touching static geometry is not evaluated inside it.
inline constexpr float kViewSurfaceOffset = 0.003f;
} // namespace placement
