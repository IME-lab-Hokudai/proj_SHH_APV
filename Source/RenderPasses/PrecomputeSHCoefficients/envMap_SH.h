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

void calculateSHCoeffs(
    std::vector<float3>& out,                // Output SH coefficients (num_basis)
    const std::vector<ProbeSampleData>& probeSamplingResults, // Probe sampling results, size = numSamples
    int numSamplePerProbe
);

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

    // this is used in testing using finite difference
void initSHBasisGradientAndHessianTables(
    const std::vector<float3>& dirSamples,
    std::vector<float>& SHBasisTableXPrime,
    std::vector<float3>& SHGradientTableXPrime,
    std::vector<float3x3>& SHHessianTableXPrime
);

void generateUniformSphereDirSamples(int sampleCount, std::vector<ProbeDirSample> &out);



float3 gradientOmega(const float3& q, const float3& n, float rInv, float cosXi, float factor);

float3x3 hessianOmega(const float3& q, const float3& n, float rInv, float cosXi, float factor);

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

void calculateGradRGBAndHessianLumSHCoeffLM(const float3& x, const ProbeSampleData* probeSamplingResults, uint32_t sampleCount, const std::vector<ProbeDirSample>& samplingDir, const int& basisIdx, GradSHCoeff& outGrad, float3x3& outHessian);
void calculateSHCoeffsGradientsRGBAndHessiansLum(std::vector<GradSHCoeff>& gradOut, std::vector<float3x3>& hessianOut, const float3& gridPos, const ProbeSampleData* probeSamplingResults,
    uint32_t sampleCount, const std::vector<ProbeDirSample>& samplingDir);
//for uniform grid
void calculateGradSHCoeffLM(const float3& x, const ProbeSampleData* probeSamplingResults, uint32_t sampleCount, const std::vector<ProbeDirSample>& samplingDir, const int& basisIdx, GradSHCoeff& outGrad);
void calculateSHCoeffsGradients(std::vector<GradSHCoeff>& gradOut, const float3& gridPos, const ProbeSampleData* probeSamplingResults, uint32_t sampleCount, const std::vector<ProbeDirSample>& samplingDir);

//progressive build related
void generateProgressiveSphereDirSamples(int sampleCount, std::vector<ProbeDirSample>& out);

void calculateSHCoeffs(
    std::vector<float3>& out,
    const ProbeSampleData* samples,
    uint32_t sampleCount
);

void calculateSHBuildMetricsOnly(
    float& coeffVecL2,
    float& maxLambdaVecL2,
    const float3& xPolar,
    const ProbeSampleData* samples,
    uint32_t sampleCount,
    const std::vector<ProbeDirSample>& samplingDirs,
    bool useRelativeMetric
);
