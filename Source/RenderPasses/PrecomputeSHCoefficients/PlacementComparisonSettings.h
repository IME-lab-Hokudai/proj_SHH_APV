#pragma once
// Settings of the placement comparison (Comparison plan: Hessian vs. adaptive probe criteria).
// Shared by PrecomputeSHCoefficients.cpp (build switch) and PlacementComparison.cpp.
#include "AdaptiveProbeVolume.h"
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

// Probe budgets = total probes of the grid, seed included (probe = grid corner, so probe count = memory).
// U8 = 9^3, U16 = 17^3, U32 = 33^3, U64 = 65^3. With the U8 seed, the 729 snapshot is the seed itself
// (a check: it must match a uniform U8 grid).
inline const std::vector<uint64_t> kProbeBudgets = { 729, 4913, 10000, 20000, 35937, 70000, 140000 };

// Leaves split per round: this fraction of the current leaf count, at most kMaxBatch.
inline constexpr float kBatchFraction = 0.01f;
inline constexpr uint32_t kMaxBatch = 256;

// Every output file starts with the scene tag: the loaded scene's file name without extension
// (e.g. SceneForData, Bistro), or kSceneTagOverride when it is not empty.
//   <scene>_BF_<score>[_S<seed>]_<budget>.txt   snapshot grid at each budget
//   <scene>_BestFirstSnapshots.csv    one row per snapshot
inline const std::string kSceneTagOverride = "";

// ---- Evaluation --------------------------------------------------------------------------------
// Per scene, in Mogwai's working directory:
//   <scene>_PlacementComparisonGrids.csv     input: one grid per line, method,type,file
//                                            (type = reference | adaptive | uniform; exactly one reference,
//                                            a uniform grid such as U64 or U128; lines starting with # skipped)
//   <scene>_PlacementComparisonSummary.csv   output: one row per grid
//   <scene>_PlacementComparisonPoints.csv    output: error per point and grid (if kWritePerPointCsv)
inline constexpr bool kWritePerPointCsv = true;

inline constexpr uint32_t kEvalResolution = 64;   // uniform lattice of cell centres per axis
inline constexpr float kRelEpsilonFraction = 0.01f; // epsilon of the percentage error, x mean reference luminance

// Run the evaluation once when the scene loads (also available as a UI button).
inline constexpr bool kEvaluateOnLoad = false;
} // namespace placement
