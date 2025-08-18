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

void initSHTable(int sh_order, int width, int height);
void decomposeSH(std::vector<float4>& out, const Falcor::ref<EnvMap>& envMap);
void reconstructSH(std::vector<float4>& sh_coeff, const Falcor::ref<EnvMap>& envMap, Falcor::ref<Device> pDevice);
float4* TranposeData(float4* data, int width, int height);

//create uniform probe grid
//TODO here is just one set of sh_coeff to test storing and loading. Need to actually find a proper way to parallax sampling
void createProbeGrid(ProbeGrid& grid, const std::vector<float4>& sh_coeff);

void saveProbeGridToFile(const ProbeGrid& grid, const std::string& path);
bool loadProbeGridFromFile( ProbeGrid& out, const std::string& path);

std::vector<ProbeDirSample> generateUniformSphereDirSamples(int sampleCount, const float3& probePos);
