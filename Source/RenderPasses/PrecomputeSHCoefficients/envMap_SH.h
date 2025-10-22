#pragma once
#include "Falcor.h"
#include "ProbeSamplingData.slang"

//using namespace std;
using namespace Falcor;

struct ProbeGrid
{
    int3 resolution;       // grid size (x, y, z)
    float3 origin;            // world-space origin of grid
    float3 spacing;           // spacing between probes
    int numBasis;
    std::vector<float4> probesSH;
    std::vector<float3> probesPos;
};

//for lat-long env map
//void initSHTable(int sh_order, int width, int height);
//void decomposeSH(std::vector<float4>& out, const Falcor::ref<EnvMap>& envMap);
//void reconstructSH(std::vector<float4>& sh_coeff, const Falcor::ref<EnvMap>& envMap, Falcor::ref<Device> pDevice);

// for probe sampling using ray tracing
void initSHTable(int sh_order, const std::vector<ProbeDirSample>& dirSamples);
void decomposeSH(
    std::vector<float4>& out,                // Output SH coefficients (num_basis)
    const std::vector<float4>& probeSamples, // Probe sampling results, size = numSamples
    int numSamplePerProbe
);
void reconstructSH(const ProbeGrid& grid, int numSamplePerProbe, std::vector<float4> & out);

// note that we will write to ylm, it is passed in as pointer to avoid repeated allocation
// calculate gradient of SH basis functions up to l = 2.
// gradient is calculated indirectly using solid spherical harmonics as described in Iwasaki sensei paper
//the code is generated using https://github.com/kiwasaki/shh_code_generator
void SHGradientL2(const float3& normDir, float* const ylm, std::array<float3, 9>& glm);
// note that we calculate up to l = 2 only 
void initSHGradientTable(const std::vector<ProbeDirSample>& dirSamples);
void computeCoeffGradients(const std::vector<ProbeDirSample>& dirSamples);

float4* TranposeData(float4* data, int width, int height);

//create uniform probe grid
void computeProbesPos(ProbeGrid& grid);

void saveProbeGridToFile(const ProbeGrid& grid, const std::string& path);
bool loadProbeGridFromFile( ProbeGrid& out, const std::string& path);

std::vector<ProbeDirSample> generateUniformSphereDirSamples(int sampleCount);
