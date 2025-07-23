#pragma once

#include <vector>
#include "Utils/Math/ScalarMath.h"
#include "Falcor.h"

using namespace std;
using namespace Falcor;

// Auxiliary function to normalize 3D vectors.
inline void normalize( float vec[ 3 ] )
{
	const float len = Falcor::math::sqrt( vec[ 0 ] * vec[ 0 ] + 
                                              vec[ 1 ] * vec[ 1 ] + 
                                              vec[ 2 ] * vec[ 2 ] );

	vec[ 0 ] /= len;
	vec[ 1 ] /= len;
	vec[ 2 ] /= len;
}
void initSHTable(int sh_order, int size);
void decomposeSH(vector<float4>& out);
void reconstructSH(std::vector<float4>& sh_coeff);

//void decomposeSH(vector<Tcolor4<float>>& out, TcubeMapTexturef& cubemap);
//void reconstructSH(vector<Tcolor4<float>>& sh_coeff, TcubeMapTexturef& cubemap);
