#include "SolidSphericalHarmonics.h"

#include "Utils/Math/MatrixTypes.h"
#include "Utils/Math/VectorMath.h"

float3 SolidSphericalHarmonics::gradOmega(float3 s, float3 x, float3 n, float N)
{
    // q = s - x
    float3 q = s - x;

    // r = ||q||
    float r = length(q);

    // cosξ = -(n · q) / r
    float cosXi = -(dot(n, q)) / r;

    // Numerical guard to avoid division by zero or near-zero values
    if (abs(cosXi) < 1e-6 || r < 1e-6)
        return float3(0.0, 0.0, 0.0);

    // Use uniform sampling with N samples
    float factor = 4.0 * M_PI / N;

    // Element-wise computation:
    // ∂xΩ_i = factor * (n_x * r + 3 * q_x * cosξ) / (r² * cosξ)
    // ∂yΩ_i = factor * (n_y * r + 3 * q_y * cosξ) / (r² * cosξ)
    // ∂zΩ_i = factor * (n_z * r + 3 * q_z * cosξ) / (r² * cosξ)
    float3 grad = factor * ((n * r + 3.0f * q * cosXi) / (r * r * cosXi));

    return grad;
}

 // Notes:
//   - gradOmega() computes ∇Ω_i (geometry term derivative).
//   - gradY_lm() provides ∇Y_l^m(ω_i) in direction space.
//   - ∂ω/∂x = -(I - ωω^T) / r maps direction derivative to spatial domain.
//   - Numerical guards in gradOmega() handle small r or cosξ cases.
//   - Typically used for SH grid refinement or Hessian computation.
//-------------------------------------------------------------------------------
float3 SolidSphericalHarmonics::gradSHCoeffLM(
    float3 x,
    const float3* s_list,
    const float3* n_list,
    const float* L_list,
    int N,
    int l,
    int m,
    const float* const SHBasisTable,
    const float3* const SHGradientTable
)
{
    float3 grad_f = float3(0.0, 0.0, 0.0);

    for (int i = 0; i < N; ++i)
    {
        float3 s = s_list[i];
        float3 n = n_list[i];
        float L = L_list[i];

        // Geometry and direction
        float3 q = s - x;
        float r = length(q);
        float3 omega_i = normalize(q); // note: omega_i points from x to s to evaluate SH basis

        // Solid angle (uniform per patch)
        float Omega_i = (4.0f * M_PI) / N;

        // Spatial derivative of Ω_i (geometry-dependent)
        float3 dOmega = gradOmega(s, x, n, N);

        // SH basis and derivative in direction space
        int numBasis = (l + 1) * (l + 1);
        int sampleIdx = i;
        int basisIdx = l * (l + 1) + m;
        float Ylm = SHBasisTable[numBasis * sampleIdx + basisIdx];
        float3 gradYlm = SHGradientTable[numBasis * sampleIdx + basisIdx];
        // Contribution of this sample to ∇f_l^m
        // Equation: ∇f_l^m ≈ Σ_i L_i [ (∇Ω_i) Y_l^m + Ω_i ∇Y_l^m ]
        grad_f += L * (dOmega * Ylm + Omega_i * gradYlm);
    }

    return grad_f;
}

float3x3 SolidSphericalHarmonics::grad2OmegaHessian(const float3& s, const float3& x, const float3& n, int N)
{
    float3 q = s - x;
    float qx = q.x, qy = q.y, qz = q.z;
    float r = length(q);

    // cosξ = -(n.q)/r
    float cosXi = -(dot(n, q)) / r;

    const float eps = 1e-6f;
    if (fabs(cosXi) < eps)
        cosXi = (cosXi < 0 ? -eps : eps);

    // Precompute powers of r
    float r2 = r * r;
    float r3 = r2 * r;
    float r4 = r2 * r2;

    // Pure second derivatives
    float d2Omega_xx = -(4.0f * M_PI / N) * (6.0f * n.x * qx * r - 3.0f * cosXi * (r2 - 5.0f * qx * qx)) / (r4 * cosXi);
    float d2Omega_yy = -(4.0f * M_PI / N) * (6.0f * n.y * qy * r - 3.0f * cosXi * (r2 - 5.0f * qy * qy)) / (r4 * cosXi);
    float d2Omega_zz = -(4.0f * M_PI / N) * (6.0f * n.z * qz * r - 3.0f * cosXi * (r2 - 5.0f * qz * qz)) / (r4 * cosXi);

    // Mixed derivatives (symmetric)
    float d2Omega_xy = -(4.0f * M_PI / N) * (((3.0f * n.x * qy + 3.0f * qx * n.y) / (r3 * cosXi)) + (15.0f * qx * qy / r4));
    float d2Omega_xz = -(4.0f * M_PI / N) * (((3.0f * n.x * qz + 3.0f * qx * n.z) / (r3 * cosXi)) + (15.0f * qx * qz / r4));
    float d2Omega_yz = -(4.0f * M_PI / N) * (((3.0f * n.y * qz + 3.0f * qy * n.z) / (r3 * cosXi)) + (15.0f * qy * qz / r4));

    // Assemble symmetric Hessian
    float3x3 H = float3x3::zeros();
    H[0][0] = d2Omega_xx;
    H[0][1] = d2Omega_xy;
    H[0][2] = d2Omega_xz;
    H[1][0] = d2Omega_xy;
    H[1][1] = d2Omega_yy;
    H[1][2] = d2Omega_yz;
    H[2][0] = d2Omega_xz;
    H[2][1] = d2Omega_yz;
    H[2][2] = d2Omega_zz;

    return H;
}

float3x3 SolidSphericalHarmonics::hessianSHCoeffLM(const float3& x,
    const float3* s_list,
    const float3* n_list,
    const float* L_list,
    int N,
    int l,
    int m,
    const float* const SHBasisTable,
    const float3* const SHGradientTable,
    const float3x3* const SHHessianTable)
{
    float3x3 H = float3x3::zeros();

    int numBasis = (l + 1) * (l + 1);
    int basisIdx = l * (l + 1) + m;
    for (int i = 0; i < N; ++i)
    {
        float3 s = s_list[i];
        float3 n = n_list[i];
        float L = L_list[i];

        float3 q = s - x;
        float r = length(q);
        float3 omega_i = normalize(q);
        float Omega_i = (4.0f * M_PI) / N;

         // Gradients and Hessian of Omega_i
        float3 dOmega = gradOmega(s, x, n, N);
        float3x3 H_Omega = grad2OmegaHessian(s, x, n, N);


    }

    return H;
}
