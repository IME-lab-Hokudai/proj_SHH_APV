#pragma once
#include "Falcor.h"

//using namespace std;
using namespace Falcor;

void initSHTable(int sh_order, int width, int height);
void decomposeSH(std::vector<float4>& out, const Falcor::ref<EnvMap>& envMap);
void reconstructSH(std::vector<float4>& sh_coeff, const Falcor::ref<EnvMap>& envMap, Falcor::ref<Device> pDevice);
float4* TranposeData(float4* data, int width, int height);

