#pragma once
#include "Falcor.h"
#include "ProbeSamplingData.slang"

//using namespace std;
using namespace Falcor;

struct GradSHCoeff
{
    float3 r; // ∇f_lm for red channel
    float3 g; // ∇f_lm for green channel
    float3 b; // ∇f_lm for blue channel
};

struct HessianSHCoeff
{
    float3x3 r; // Hessian of f_lm for red channel
    float3x3 g; // Hessian of f_lm for green channel
    float3x3 b; // Hessian of f_lm for blue channel
};

struct ProbeGrid
{
    int3 resolution;       // grid size (x, y, z)
    float3 origin;            // world-space origin of grid
    float3 spacing;           // spacing between probes
    int numBasis;
    std::vector<float4> probesSHCoeffs;
    std::vector<GradSHCoeff> probesSHCoeffsGradients;
    std::vector<HessianSHCoeff> probesSHCoeffsHessians;
    std::vector<float3> probesPos;
};

//for lat-long env map
//void initSHTable(int sh_order, int width, int height);
//void decomposeSH(std::vector<float4>& out, const Falcor::ref<EnvMap>& envMap);
//void reconstructSH(std::vector<float4>& sh_coeff, const Falcor::ref<EnvMap>& envMap, Falcor::ref<Device> pDevice);

// for probe sampling using ray tracing
void initSHTable(int sh_order, const std::vector<ProbeDirSample>& dirSamples);
void calculateSHCoeffs(
    std::vector<float3>& out,                // Output SH coefficients (num_basis)
    const std::vector<ProbeSampleData>& probeSamplingResults, // Probe sampling results, size = numSamples
    int numSamplePerProbe
);
void calculateSHCoeffsGradientsAndHessians(
    std::vector<GradSHCoeff>& gradOut,
    std::vector<HessianSHCoeff>& hessianOut,
    const float3& gridPos,
    const std::vector<ProbeSampleData>& probeSamplingResults,
    const std::vector<ProbeDirSample>& samplingDir
);



void reconstructSH(const ProbeGrid& grid, int numSamplePerProbe, std::vector<float4> & out);

// calculate gradient and hessian of SH basis functions up to l = 2.
// gradient is calculated indirectly using solid spherical harmonics as described in Iwasaki sensei paper
//the code is generated using https://github.com/kiwasaki/shh_code_generator
// note: it is for SH basis (Ylm) NOT for SH coefficients (Flm). Hlm and Glm are used in equation 2 and 4
// it computes Glm and Hlm wrt direction vector (not spatial position)
void SHGradientAndHessianL2(const float3& normDir, std::array<float, 9>& ylm, std::array<float3, 9>& glm, std::array<float3x3, 9>& hlm);
void SHGradientL2(const float3& normDir, std::array<float, 9>& ylm, std::array<float3, 9>& glm);

// note that we calculate up to l = 2 only
void initSHBasisGradientAndHessianTables(const std::vector<ProbeDirSample>& dirSamples);

/** * Converts linearized basis index to (l, m).
 * Assumes encoding: basisIdx = l * (l + 1) + m
 * This works because the index range for band l is [l^2, (l+1)^2 - 1].
 */
void getLMFromBasisIdx(int basisIdx, int& l, int& m);

float3 computeKrivanekBasisGradient(int l, int m, float r, float3 dir, float evalYlmminus);
void computeKrivanekCoeffLMGradient(float3 x, std::vector<ProbeSampleData> samplingData, int basisIdx, float3& outGrad);
    // this is used in testing using finite difference
void initSHBasisGradientAndHessianTables(
    const std::vector<float3>& dirSamples,
    std::vector<float>& SHBasisTableXPrime,
    std::vector<float3>& SHGradientTableXPrime,
    std::vector<float3x3>& SHHessianTableXPrime
);

float4* TranposeData(float4* data, int width, int height);

//create uniform probe grid
void computeProbesPos(ProbeGrid& grid);

void saveProbeGridToFile(const ProbeGrid& grid, const std::string& path);
void saveProbeGridToFileWithGradAndHessian(const ProbeGrid& grid, const std::string& path);
bool loadProbeGridFromFileWithGradAndHessian(ProbeGrid& out, const std::string& path);
bool loadProbeGridFromFile( ProbeGrid& out, const std::string& path);

void generateUniformSphereDirSamples(int sampleCount, std::vector<ProbeDirSample> &out);

//-------------------------------------------------------------------------------
// Function: gradOmega
// Purpose : Compute the spatial gradient of Ω_i with respect to the grid point x.
//
// Equation 3:
//   ∂Ω_i/∂x = (4π / N) * ((n_x * r + 3 * q_x * cosξ) / (r² * cosξ))
//   (and similarly for ∂Ω_i/∂y, ∂Ω_i/∂z)
//
// Input Parameters:
//   s  : float3 - Position of the environmental patch (source point on the scene surface)
//   x  : float3 - Position of the grid point where gradient of Ω_i is evaluated
//   n  : float3 - Surface normal at the environmental patch
//   N  : int  - number of samples
//
// Intermediate Terms:
//   q      = s - x                     // Vector from grid point to patch
//   r      = ||q||                     // Distance between grid point and patch
//   cosξ   = -(n · q) / r              // Cosine of angle between n and q
//
// Output:
//   Returns float3(∂Ω_i/∂x, ∂Ω_i/∂y, ∂Ω_i/∂z)
//   - The spatial gradient of Ω_i evaluated at grid point x.
//
// Notes:
//   - Used for precomputed SH grid construction or adaptive refinement.
//   - Avoids division by near-zero cosξ or r using numerical guards.
//   - The direction of q (from x → s) follows the convention of incoming radiance.
//-------------------------------------------------------------------------------

float3 gradientOmega(float3 s, float3 x, float3 n, float N);


//-------------------------------------------------------------------------------
// Compute the full 3x3 Hessian of the solid angle Ω_i
//-------------------------------------------------------------------------------
// Inputs:
//   s  : float3 - Position of the environmental patch (source point on the scene surface)
//   x  : float3 - Position of the grid point where gradient of Ω_i is evaluated
//   n  : float3 - Surface normal at the environmental patch
//   N  : int  - number of samples
// Equation 5 and 6
// ∂_xx Ω_i= -∂_(q_x ) ∂_(q_x ) Ω_i= -4π/N⋅(6n_x q_x r-3cosξ(r^2- 5q_x^2))/(r^4 cosξ)
// ∂_yx Ω_i = -∂_(q_y) ∂_(q_x) Ω_i = -4π / N((3n_x q_y + 3q_x n_y) / (r ^ 3 cosξ) + (15q_x q_y) / r ^ 4)
// Returns:
//    3x3 Hessian matrix H, where H(i,j) = ∂²Ω_i / ∂x_i ∂x_j
//
// Notes:
//   - q = s - x
//   - r = length(q)
//   - cosξ = -(n ⋅ q)/r
//   - Hessian contains both pure derivatives (∂²/∂x², ∂²/∂y², ∂²/∂z²)
//     and mixed derivatives (∂²/∂x∂y, etc.)
//   - Be careful of singularities when cosξ → 0; clamp with a small epsilon
//-------------------------------------------------------------------------------
float3x3 hessianOmega(const float3& s, const float3& x, const float3& n, int N);

//-------------------------------------------------------------------------------
// Function: grad_SH
// Purpose : Compute the full spatial gradient ∇f_l^m at a grid point.
//
// Equation 2:
//   ∇f_l^m ≈ Σ_i L(ω_i) [ (∇Ω_i) * Y_l^m(ω_i) + Ω_i * ∇Y_l^m(ω_i) ]
//
// Input Parameters:
//   s_list : pointer to float3 array - Positions of environmental patches
//   n_list : pointer to float3 array - Normals of environmental patches
//   L_list : pointer to float array  - Radiance (or irradiance) per patch
//   x      : float3                  - Grid point position
//   N  : int  - number of samples
//   l, m   : int                     - SH band and order indices
//
// Intermediate Terms:
//   q      = s - x                   // Vector from grid point to patch
//   r      = ||q||                   // Distance between grid point and patch
//   ω_i    = normalize(q)            // Direction from grid to patch
//   cosξ   = -(n · q) / r            // Cosine of incident angle
//
// Output:
//   Returns float3(∂_x f_l^m, ∂_y f_l^m, ∂_z f_l^m)
//   - The spatial gradient of f_l^m evaluated at grid point x.
void calculateGradAndHessianSHCoeffLM(
    const float3 &x,
    const std::vector<ProbeSampleData>& samplingData,
    const std::vector<ProbeDirSample>& samplingDir,
    const int &basisIdx,
    GradSHCoeff& outGrad,
    HessianSHCoeff& outHessian
);

//for uniform grid
void calculateGradSHCoeffLM(const float3& x, const std::vector<ProbeSampleData>& samplingData, const std::vector<ProbeDirSample>& samplingDir, const int& basisIdx, GradSHCoeff& outGrad);
void calculateSHCoeffsGradients(std::vector<GradSHCoeff>& gradOut, const float3& gridPos, const std::vector<ProbeSampleData>& probeSamplingResults, const std::vector<ProbeDirSample>& samplingDir);
//-------------------------------------------------------------------------------
// Compute full 3x3 Hessian of SH coefficient f_l^m at a grid point
//-------------------------------------------------------------------------------
// Inputs:
//   x                 : float3, grid point
//   s_list            : array of float3, patch/sample positions
//   n_list            : array of float3, patch normals
//   L_list            : array of float, patch radiance/intensity
//   N                 : int, number of samples
//   l, m              : SH band and order
//   SHBasisTable      : float*, precomputed Y_l^m for all samples
//   SHGradientTable   : float3*, precomputed ∇Y_l^m for all samples
//   SHHessianTable    : glm::mat3*, precomputed ∇²Y_l^m for all samples
//
// Returns:
//   full Hessian ∂² f_l^m / ∂x_j ∂x_k
//-------------------------------------------------------------------------------
HessianSHCoeff hessianSHCoeffLM(
    const float3& x,
    const std::vector<ProbeSampleData>& samplingData,
    const int& basisIdx
);


// ****functions for verifying SH gradient and hessian calculation using numerical finite differentiation****

// Generate a single vector containing all positions needed for ∂f/∂x finite difference
// Order per sample: [center, x+h, x-h]
std::vector<float3> generateVerificationPositions(float y = 0.2f, float extent = 0.25f, uint32_t resolution = 32, float h = 0.02f);
// Order per sample: [center, x+h z+h, x+h z-h, x-h z+k, x-h z-k]
std::vector<float3> generateVerificationPositionsMixed(float y = 0.2f, float extent = 0.25f, uint32_t resolution = 32, float h = 0.02f);

// we evaluate coeff for just one basis and channel R for verification purpose
float calculateChannelRSHCoeffLM(int basisIdx, const std::vector<ProbeSampleData>& probeSamplingResults);
void calculateChannelRGradAndHessianSHCoeffLM(
    const float3& x,
    const std::vector<ProbeSampleData>& samplingData,
    const std::vector<ProbeDirSample>& samplingDir,
    const int& basisIdx,
    float3& outGrad,
    float3x3& outHessian
);

float3 testComputeSHGrad();

void calculateGradRGBAndHessianLumSHCoeffLM(const float3& x, const std::vector<ProbeSampleData>& samplingData, const std::vector<ProbeDirSample>& samplingDir, const int& basisIdx, GradSHCoeff& outGrad, float3x3& outHessian);

void calculateSHCoeffsGradientsRGBAndHessiansLum(std::vector<GradSHCoeff>& gradOut, std::vector<float3x3>& hessianOut, const float3& gridPos, const std::vector<ProbeSampleData>& probeSamplingResults, const std::vector<ProbeDirSample>& samplingDir);

void calculateSHCoeffsRadialMoments(
    std::vector<float>& outMean,      // Output: 9 coeffs for E[r]
    std::vector<float>& outMeanSq,    // Output: 9 coeffs for E[r^2]
    const std::vector<ProbeSampleData>& probeSamplingResults,
    int numSamplePerProbe
);
