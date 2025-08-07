#pragma once
#include "Falcor.h"

//using namespace std;
using namespace Falcor;

struct SHProbe
{
    std::vector<float4> shCoeffs; // RGBA SH coefficients (L2)
};

struct ProbeGrid
{
    int3 resolution;       // grid size (x, y, z)
    float pad0 = 0.0f; //padding
    float3 origin;            // world-space origin of grid
    float pad1 = 0.0f;
    float3 spacing;           // spacing between probes
    float pad2 = 0.0f;
    std::vector<SHProbe> probes;
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

