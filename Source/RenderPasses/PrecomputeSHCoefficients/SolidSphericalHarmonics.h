#pragma once

#include <vector>

#include "Utils/Math/MatrixTypes.h"
#include "Utils/Math/VectorTypes.h"
using namespace std;
using namespace Falcor;

#ifndef M_PI
#define M_PI 3.14159268
#endif

class SolidSphericalHarmonics
{
public:
private:

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

    float3 gradOmega(float3 s, float3 x, float3 n, float N);

//-------------------------------------------------------------------------------
    // Function: grad_SH
    // Purpose : Compute the full spatial gradient ∇f_l^m at a grid point.
    //
    // Equation:
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
    //

    float3 gradSHCoeffLM(
        float3 x,
        const float3* s_list,
        const float3* n_list,
        const float* L_list,
        int N,
        int l,
        int m,
        const float* const SHBasisTable,
        const float3* const SHGradientTable
    );

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
    //∂_yx Ω_i = -∂_(q_y) ∂_(q_x) Ω_i = -4π / N((3n_x q_y + 3q_x n_y) / (r ^ 3 cosξ) + (15q_x q_y) / r ^ 4)
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
    float3x3 grad2OmegaHessian(const float3& s, const float3& x, const float3& n, int N);

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
    float3x3 hessianSHCoeffLM(
        const float3& x,
        const float3* s_list,
        const float3* n_list,
        const float* L_list,
        int N,
        int l,
        int m,
        const float* const SHBasisTable,
        const float3* const SHGradientTable,
        const float3x3* const SHHessianTable
    );
};
