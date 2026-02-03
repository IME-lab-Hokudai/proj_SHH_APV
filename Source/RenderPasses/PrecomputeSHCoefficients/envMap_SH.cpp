#include "sphericalHarmonics.h"
#include "CommonDefines.h"
#include "envMap_SH.h"

#include <fstream>
#include <random>

float*  dOmega;
float*  SHBasisTable;
float3* SHGradientTable;
float3x3* SHHessianTable;
int shOrder = -1;

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "stb_image.h"
#include "stb_image_write.h"
//-----------------------------------------------------------------------
//void sphericalHarmonics16( const float x, const float y, const float z, float sh[ 16 ] )
//{
//	// A is normalization factor: A = sqrt( ( 2 * l + 1 ) / ( 4 * PI ) ) )
//	float c[ 16 ];
//
//	c[ 0] =  0.282095f;
//	c[ 1] = -0.488603f;
//	c[ 2] =  0.488603f;
//	c[ 3] = -0.488603f;
//
//	c[ 4] =  1.092548f;
//	c[ 5] = -1.092548f;
//	c[ 6] =  0.315392f;
//	c[ 7] = -1.092548f;
//	c[ 8] =  0.546274f;
//
//	c[ 9] = -0.590044f;			// Y{3, -3}
//	c[10] =  2.89061f;			// Y{3, -2}
//	c[11] = -0.457046f;			// Y{3, -1}
//	c[12] =  0.373176f;			// Y{3,  0}
//	c[13] = -0.457046f;			// Y{3,  1}
//	c[14] =  1.44531f;			// Y{3,  2}
//	c[15] = -0.590044f;			// Y{3,  3}
//
//	//for(int i=0; i<16; i++) c[i] = 1.0;
//
//	// Y_{0,0} 
//	sh[ 0] = c[ 0];
//
//	// Y_{1,-1} Y_{1,0}, Y_{1,1}
//	sh[ 1] = c[ 1] * y;
//	sh[ 2] = c[ 2] * z;
//	sh[ 3] = c[ 3] * x;
//
//	// Y_{2, -2} Y_{2,-1}, Y_{2,1}
//	sh[ 4] = c[ 4] * x * y;
//	sh[ 5] = c[ 5] * y * z;
//	sh[ 7] = c[ 7] * x * z;
//
//	// Y_{2,0} 
//	sh[ 6] = c[ 6] * (3.0f * z * z - 1.0f);
//
//	// Y_{2,2} 
//	sh[ 8] = c[ 8] * (x * x - y * y);
//
//	// Y_{3, -3} = A * sqrt(5/8) * (3 * x^2 * y - y^3)
//	sh[ 9] = c[ 9] * (3.0f * x * x * y - y * y * y); 
//
//	// Y_{3, -2} = A * sqrt(15) * x * y * z 
//	sh[10] = c[10] * x * y * z;
//
//	// Y_{3, -1} = A * sqrt(3/8) * y * (5 * z^2 - 1)
//	sh[11] = c[11] * y * (5.0f * z * z - 1.0f);
//
//	// Y_{3,  0} = A * (1/2) * (5 * z^3 - 3 *z)	
//	sh[12] = c[12] * (5.0f * z * z * z - 3 * z);
//
//	// Y_{3,  1} = A * sqrt(3/8) * x * (5 * z^2 - 1)
//	sh[13] = c[13] * x * (5.0f * z * z  - 1.0f);
//
//	// Y_{3,  2} = A * sqrt(15/4) * z *(x^2 - y^2)
//	sh[14] = c[14] * z * (x * x - y * y);
//
//	// Y_{3,  3} = A * sqrt(5/8) * (x^3 - 3 * x * y^2)
//	sh[15] = c[15] * (x * x * x - 3.0f * x * y * y);
//}

//  precomputed pixel solid angle weight sin(theta)*dtheta*dphi
//void calcDeltaFormFactorEquirect(int width, int height)
//{
//    dOmega = new float[width * height];
//    double dTheta = M_PI / height;
//    double dPhi = 2.0 * M_PI / width;
//
//    for (int y = 0; y < height; ++y)
//    {
//        // theta from 0 to pi
//        double theta = M_PI * (y + 0.5) / height;
//        double sinTheta = sin(theta);
//
//        for (int x = 0; x < width; ++x)
//        {
//            dOmega[TWO_D_TO_ONE_D(x,y,width)] = (float)sinTheta * dTheta * dPhi;
//        }
//    }
//
//    //double sum = 0.0;
//    //for (int y = 0; y < height; ++y)
//    //    for (int x = 0; x < width; ++x)
//    //        sum += dOmega[y * width + x];
//
//    //std::cout << "Total Omega = " << sum << std::endl; // Should be ≈ 4π
//}

//this function precalculate value of all SH basis in all direction and store in SHTableLookup 
//then when we precompute coefficients we look up from SHTableLookup to calculate Coeff = dOmega * env_map * Y
//void initSHTable(int sh_order, int width, int height)
//{
//    calcDeltaFormFactorEquirect(width, height);
//    shOrder = sh_order;
//    SphericalHarmonics sh(sh_order);
//    int num_basis = sh.getNumBasis();
//    // Allocate SH basis table for each pixel (flattened 3D array: [basis][width][height])
//    SHBasisTable = new float[num_basis * width * height];
//
//    // Preallocate vector for SH basis
//    vector<double> y(num_basis);
//    for (int y_idx = 0; y_idx < height; ++y_idx)
//    {
//        // θ = latitude angle from 0 (top) to π (bottom)
//        double v = (y_idx + 0.5) / height;
//        double theta = v * M_PI;
//
//        for (int x_idx = 0; x_idx < width; ++x_idx)
//        {
//            // φ = longitude angle from 0 to 2π
//            double u = (x_idx + 0.5) / width;
//            double phi = u * 2.0 * M_PI;
//
//            // Compute SH basis at (θ, φ)
//            sh.calcSHBasis(y, cos(theta), phi); // ct = cos(θ)
//
//            //storing all Ylm for this direction (i.e 9 Y if l = 2)
//            for (int i = 0; i < num_basis; ++i)
//            {
//                SHBasisTable[THREE_D_TO_ONE_D(i, x_idx, y_idx, width , height)] = y[i];
//            }
//        }
//    }
//    //float Y00 = envSHTable[THREE_D_TO_ONE_D(0, width / 2, height / 2, width, height)];
//    //std::cout << "Y00 = " << Y00 << std::endl;
//}

//data layout
//sample0: b0 b1 b2 ... bN
//sample1 : b0 b1 b2... bN
//sample2 : b0 b1 b2... bN
//--> indexing scheme SHBasisTable[num_basis * sampleIdx + basisIdx]
void initSHTable(int sh_order, const std::vector<ProbeDirSample>& dirSamples)
{
    shOrder = sh_order;
    SphericalHarmonics sh(sh_order);
    int num_basis = sh.getNumBasis();
    int numSamplePerProbe = (int)dirSamples.size();

    // Allocate SH basis table for each direction (flattened 2D array: [basis][sampleIdx])
    if (SHBasisTable)
        delete[] SHBasisTable;
    SHBasisTable = new float[num_basis * numSamplePerProbe];

    // Preallocate vector for SH basis
    std::vector<double> y(num_basis);

    for (int sampleIdx = 0; sampleIdx < numSamplePerProbe; ++sampleIdx)
    {
        const float3& dir = dirSamples[sampleIdx].dir;
        //const float3& dir = -dirSamples[sampleIdx].dir;

        // Convert to spherical coordinates
        double theta = acos(clamp(dir.y, -1.0f, 1.0f)); // y_ is cos(theta)
        double phi = atan2(dir.z, dir.x);                         // Note: atan2(z, x) for Falcor's Y-up

        // Compute SH basis for this direction
        sh.calcSHBasis(y, cos(theta), phi);

        for (int basisIdx = 0; basisIdx < num_basis; ++basisIdx)
        {
            //SHBasisTable[basisIdx * numSamplePerProbe + sampleIdx] = (float)y[basisIdx];
            SHBasisTable[num_basis * sampleIdx + basisIdx] = (float)y[basisIdx];
        }
    }
}

void initSHBasisGradientAndHessianTables(const std::vector<ProbeDirSample>& dirSamples)
{
    shOrder = 2; // fixed to L2 for now
    int numBasis = 9;
    std::array<float3, 9> glm;
    std::array<float3x3, 9> hlm;
    std::array<float, 9> ylm; // local storage for SH values
    int numSamplesPerProbe = (int)dirSamples.size();

    // Allocate SH basis table for each direction (flattened 2D array: [basis][sampleIdx])
    if (SHBasisTable)
        delete[] SHBasisTable;
    SHBasisTable = new float[9 * numSamplesPerProbe];

    if (SHGradientTable)
        delete[] SHGradientTable;
    SHGradientTable = new float3[9 * numSamplesPerProbe];

      if (SHHessianTable)
        delete[] SHHessianTable;
      SHHessianTable = new float3x3[9 * numSamplesPerProbe];

    for (int sampleIdx = 0; sampleIdx < numSamplesPerProbe; ++sampleIdx)
    {
        const float3& dir = dirSamples[sampleIdx].dir;

        // compute SH gradients for this direction
        SHGradientAndHessianL2(dir, ylm, glm, hlm);

        // store into SHGradientTable basisidx = l(l+1)+m
        for (int basisIdx = 0; basisIdx < numBasis; ++basisIdx)
        {
            SHBasisTable[sampleIdx * numBasis + basisIdx] = ylm[basisIdx];
            SHGradientTable[sampleIdx * numBasis + basisIdx] = glm[basisIdx];
            SHHessianTable[sampleIdx * numBasis + basisIdx] = hlm[basisIdx];
        }
    }
}

void getLMFromBasisIdx(int basisIdx, int& l, int& m)
{
    // l = floor(sqrt(basisIdx))
    l = static_cast<int>(std::sqrt(static_cast<float>(basisIdx)));
    // m = basisIdx - CenterOfBand(l)
    m = basisIdx - (l * (l + 1));
}

float3 computeKrivanekBasisGradient(int l, int m, float r, float3 dir, float evalYlmminus)
{
    // --- Step A: Spherical Coordinates ---
    float z_safe = std::clamp(dir.z, -1.0f, 1.0f);
    double theta = std::acos(z_safe);
    double phi = std::atan2(dir.y, dir.x);
    if (phi < 0.0)
        phi += 2.0 * M_PI;
    // x = sintheta*cosphi | y = sintheta*sinphi | z = costheta (r = 1 so omitted)
    double cos_theta = z_safe;
    double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
    double cos_phi = (sin_theta > 1e-6) ? (dir.x / sin_theta) : 1.0;
    double sin_phi = (sin_theta > 1e-6) ? (dir.y / sin_theta) : 0.0;

        // 1. dYlm/dTheta
     double dYlm_dTheta = 0.0;
    
    double K = getNormalizationK(l, m);
    double dP_dx_val = get_dP_dx(l, m, cos_theta);

    if (m > 0)
        dYlm_dTheta = -std::sqrt(2.0) * K * std::cos(m * phi) * sin_theta * dP_dx_val;
    else if (m < 0)
        dYlm_dTheta = -std::sqrt(2.0) * K * std::sin(-m * phi) * sin_theta * dP_dx_val;
    else
        dYlm_dTheta = -K * sin_theta * dP_dx_val;
    

    // 2. dY/dPhi ( = -m * Y_l^{-m} )
    double dYlm_dPhi = (m == 0) ? 0.0 : -(double)m * evalYlmminus;

    // --- Step C: Spatial Gradients (Chain Rule) ---
    const float EPSILON = 1e-6f;
    float r_safe = (r < EPSILON) ? EPSILON : r;
    float inv_r = 1.0f / r_safe;
    float inv_r_sin_theta = (sin_theta > EPSILON) ? (inv_r / (float)sin_theta) : 0.0f;

    // Geometric Derivatives (Eq 9 & 11)
    float dTheta_dx = -(float)(cos_theta * cos_phi) * inv_r;
    float dTheta_dy = -(float)(cos_theta * sin_phi) * inv_r;
    float dPhi_dx = (float)sin_phi * inv_r_sin_theta;
    float dPhi_dy = -(float)cos_phi * inv_r_sin_theta;

    // Final Chain Rule (Eq 8)
    float3 grad;
    grad.x = (float)(dTheta_dx * dYlm_dTheta + dPhi_dx * dYlm_dPhi);
    grad.y = (float)(dTheta_dy * dYlm_dTheta + dPhi_dy * dYlm_dPhi);
    grad.z = 0.0f;

    return grad;
}

void computeKrivanekCoeffLMGradient(float3 x, std::vector<ProbeSampleData> samplingData, int basisIdx, float3& outGrad)
{
    outGrad = float3(0.0f, 0.0f, 0.0f);
    int samplingSize = samplingData.size();
    int numBasis = 9;
    int l, m;
    getLMFromBasisIdx(basisIdx, l, m);
    int basisIdxMMinus = l * (l + 1) + (-m);
    for (int sampleIdx = 0; sampleIdx < samplingSize; ++sampleIdx)
    {
        if (samplingData[sampleIdx].hitT < 0.0f) // ray miss
            continue;
        float3 s = float3(samplingData[sampleIdx].s.x, samplingData[sampleIdx].s.y, samplingData[sampleIdx].s.z);
        float3 n = float3(samplingData[sampleIdx].n.x, samplingData[sampleIdx].n.y, samplingData[sampleIdx].n.z);
        float3 L = float3(samplingData[sampleIdx].Li.x, samplingData[sampleIdx].Li.y, samplingData[sampleIdx].Li.z);

        // Solid angle (uniform per patch)
        float Omega_i = (4.0f * M_PI) / (float)samplingSize;
        // Gradient ∂_x Ω_i
        float3 gradOmega = gradientOmega(s, x, n, samplingSize);
        // q = s - x
        float3 q = s - x;

        // r = ||q||
        float r = length(q);
        

        float Ylm = SHBasisTable[sampleIdx * numBasis + basisIdx];
        float Ylmminus = SHBasisTable[sampleIdx * numBasis + basisIdxMMinus];

        float3 gradYlm = computeKrivanekBasisGradient(l, m, r, normalize(q), Ylmminus);

        float3 contrib = gradOmega * Ylm + Omega_i * gradYlm;
        outGrad += (L.r * contrib);
    }
}

void initSHBasisGradientAndHessianTables(const std::vector<float3>& dirSamples,
                                         std::vector<float>&  SHBasisTableXPrime,
                                         std::vector<float3>&  SHGradientTableXPrime,
                                         std::vector<float3x3>&  SHHessianTableXPrime)
{
    int numBasis = 9;
    std::array<float3, 9> glm;
    std::array<float3x3, 9> hlm;
    std::array<float, 9> ylm; // local storage for SH values
    int numSamplesPerProbe = (int)dirSamples.size();

    for (int sampleIdx = 0; sampleIdx < numSamplesPerProbe; ++sampleIdx)
    {
        const float3& dir = dirSamples[sampleIdx];

        // compute SH gradients for this direction
        SHGradientAndHessianL2(dir, ylm, glm, hlm);

        // store into SHGradientTable
        for (int basisIdx = 0; basisIdx < numBasis; ++basisIdx)
        {
            SHBasisTableXPrime[sampleIdx * numBasis + basisIdx] = ylm[basisIdx];
            SHGradientTableXPrime[sampleIdx * numBasis + basisIdx] = glm[basisIdx];
            SHHessianTableXPrime[sampleIdx * numBasis + basisIdx] = hlm[basisIdx];
        }
    }
}

void calculateSHCoeffs(
    std::vector<float3>& out,                // Output SH coefficients (num_basis)
    const std::vector<ProbeSampleData>& probeSamplingResults, // Probe sampling results, size = numSamples
    int numSamplePerProbe
)
{
    if (shOrder == -1 || SHBasisTable == nullptr)
    {
        logError("call initSHTable before calculate SH coeffs!");
        return;
    }

    int num_basis = 9;
    out.clear();
    out.resize(num_basis);

    float weight = 4.0f * float(M_PI) / float(numSamplePerProbe); //because uniform sampling

    // For each SH basis
    for (int basisIdx = 0; basisIdx < num_basis; ++basisIdx)
    {
        double r = 0.0, g = 0.0, b = 0.0, a = 0.0;

        // For each direction sample
        for (int sampleIdx = 0; sampleIdx < numSamplePerProbe; ++sampleIdx)
        {
            const float4& sample = probeSamplingResults[sampleIdx].Li;
            float shBasis = SHBasisTable[num_basis * sampleIdx + basisIdx]; // SH basis value for this direction

            r += (double)sample.x * shBasis * weight;
            g += (double)sample.y * shBasis * weight;
            b += (double)sample.z * shBasis * weight;
            //a += (double)sample.w * shBasis * weight;
        }
        //out[basisIdx] = float3((float)r, (float)g, (float)b, (float)a);
        out[basisIdx] = float3((float)r, (float)g, (float)b);
    }
}

void calculateSHCoeffsGradientsAndHessians(
    std::vector<GradSHCoeff>& gradOut,
    std::vector<HessianSHCoeff>& hessianOut,
    const float3& gridPos,
    const std::vector<ProbeSampleData>& probeSamplingResults,
    const std::vector<ProbeDirSample>& samplingDir
)
{
    gradOut.clear();
    gradOut.resize(9);
    hessianOut.clear();
    hessianOut.resize(9);
    // For each SH basis
    for (int basisIdx = 0; basisIdx < 9; ++basisIdx)
    {
        GradSHCoeff grad;
        HessianSHCoeff hessian;
        calculateGradAndHessianSHCoeffLM(gridPos, probeSamplingResults, samplingDir, basisIdx, grad, hessian);
        gradOut[basisIdx] = grad;
        hessianOut[basisIdx] = hessian;
    }
}

void calculateSHCoeffsGradients(
    std::vector<GradSHCoeff>& gradOut,
    const float3& gridPos,
    const std::vector<ProbeSampleData>& probeSamplingResults,
    const std::vector<ProbeDirSample>& samplingDir
)
{
    gradOut.clear();
    gradOut.resize(9);
    // For each SH basis
    for (int basisIdx = 0; basisIdx < 9; ++basisIdx)
    {
        GradSHCoeff grad;
        HessianSHCoeff hessian;
        calculateGradSHCoeffLM(gridPos, probeSamplingResults, samplingDir, basisIdx, grad);
        gradOut[basisIdx] = grad;
    }
}

// irr = sumlm(coefflm* Ylm)
void reconstructSH(const ProbeGrid& grid, int numSamplePerProbe, std::vector<float4>& reconstructedData)
{
    int num_basis = grid.numBasis;
    int numProbes = grid.resolution.x * grid.resolution.y * grid.resolution.z;

    reconstructedData.clear();
    reconstructedData.resize(numProbes * numSamplePerProbe);

    for (int probeIdx = 0; probeIdx < numProbes; ++probeIdx)
    {
        for (int sampleIdx = 0; sampleIdx < numSamplePerProbe; ++sampleIdx)
        {
            float4 irr = float4(0, 0, 0, 0);
            for (int basisIdx = 0; basisIdx < num_basis; ++basisIdx)
            {
                irr += grid.probesSHCoeffs[probeIdx * num_basis + basisIdx] * SHBasisTable[num_basis * sampleIdx + basisIdx];
            }
            reconstructedData[probeIdx * numSamplePerProbe + sampleIdx] = irr;
        }
    }
}

void SHGradientAndHessianL2(const float3& normDir, std::array<float, 9>& ylm, std::array<float3, 9>& glm, std::array<float3x3, 9>& hlm)
{
    float x = normDir.x;
    float y = normDir.y;
    float z = normDir.z;
    float c0, c1, cm, cs, s0, s1, sm, ss, tmp, tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, lx, ly, lz;
    const float x2 = x * x;
    const float y2 = y * y;
    const float z2 = z * z;
    const float xy = x * y;
    const float xz = x * z;
    const float yz = y * z;
    std::array<float, 9> qlm{};
    // zonal harmonics(m=0)
    ylm[0] = qlm[0] = 0.2820947917738781f;
    ylm[2] = qlm[3] = 0.4886025119029199f * z;
    ylm[6] = qlm[6] = 0.9461746957575601f * z2 - 0.3153915652525200f;
    c0 = x;
    s0 = y;
    c1 = 1;
    s1 = 0;
    // m = 001
    qlm[4] = -0.4886025119029200f;
    ylm[3] = qlm[4] * c0;
    ylm[1] = qlm[4] * s0;
    qlm[7] = -1.0925484305920792f * z;
    ylm[7] = qlm[7] * c0;
    ylm[5] = qlm[7] * s0;
    s1 = s0;
    c1 = c0;
    c0 = x * c1 - y * s1;
    s0 = y * c1 + x * s1;
    // m = 002
    qlm[8] = 0.5462742152960396f;
    ylm[8] = qlm[8] * c0;
    ylm[4] = qlm[8] * s0;
    // calculate gradient
    glm[0] = {};
    hlm[0] =float3x3::zeros();
    lx = x;
    ly = y;
    lz = z;
    glm[2].x = -x * ylm[2];
    glm[2].y = -y * ylm[2];
    glm[2].z = -z * ylm[2] + 0.4886025119029199f;
    hlm[2][0][0] = -((1 - x2) * ylm[2] + 2 * x * glm[2].x);
    hlm[2][0][1] = -(-xy * ylm[2] + y * glm[2].x + x * glm[2].y);
    hlm[2][0][2] = -(-xz * ylm[2] + z * glm[2].x + x * glm[2].z);
    hlm[2][1][0] = hlm[2][0][1];
    hlm[2][1][1] = -((1 - y2) * ylm[2] + 2 * y * glm[2].y);
    hlm[2][1][2] = -(-yz * ylm[2] + z * glm[2].y + y * glm[2].z);
    hlm[2][2][0] = hlm[2][0][2];
    hlm[2][2][1] = hlm[2][1][2];
    hlm[2][2][2] = -((1 - z2) * ylm[2] + 2 * z * glm[2].z);
    lx += x;
    ly += y;
    lz += z;
    tmp = 2 * ylm[6];
    tmp0 = 1.2909944487358056f * qlm[4] - tmp;
    glm[6].x = x * tmp0;
    glm[6].y = y * tmp0;
    glm[6].z = -z * tmp + 2.5819888974716116f * qlm[3];
    tmp1 = 0 * tmp;
    tmp2 = 1.5811388300841898f * qlm[2] - tmp1;
    tmp3 = 1.5811388300841898f * qlm[1] - glm[6].z;
    hlm[6][0][0] = x2 * tmp2 + tmp0 - 2 * lx * glm[6].x;
    hlm[6][0][1] = hlm[6][1][0] = xy * tmp2 - ly * glm[6].x - lx * glm[6].y;
    hlm[6][1][1] = y2 * tmp2 + tmp0 - 2 * ly * glm[6].y;
    hlm[6][0][2] = hlm[6][2][0] = -xz * tmp1 + lx * tmp3 - lz * glm[6].x;
    hlm[6][1][2] = hlm[6][2][1] = -yz * tmp1 + ly * tmp3 - lz * glm[6].y;
    hlm[6][2][2] = -z2 * tmp1 - 2 * lz * glm[6].z + 4.4721359549995796f * qlm[0] - tmp;
    lx += x;
    ly += y;
    lz += z;
    c0 = x;
    s0 = y;
    c1 = 1;
    s1 = 0;
    cm = 1;
    sm = 0;
    cs = 0;
    ss = 0;
    // m = 001
    lx = 1 * x;
    ly = 1 * y;
    lz = 1 * z;
    tmp0 = qlm[4] * cm;
    tmp1 = qlm[4] * sm;
    glm[3].x = -lx * ylm[3] + tmp0;
    glm[3].y = -ly * ylm[3] - tmp1;
    glm[3].z = -lz * ylm[3];
    glm[1].x = -lx * ylm[1] + tmp1;
    glm[1].y = -ly * ylm[1] + tmp0;
    glm[1].z = -lz * ylm[1];
    tmp0 = 1 * ylm[3];
    tmp1 = -1 * tmp0;
    tmp2 = qlm[4] * cs;
    tmp3 = qlm[4] * ss;
    hlm[3][0][0] = -tmp1 * x2 - tmp0 + tmp2 - (lx + lx) * glm[3].x;
    hlm[3][0][1] = hlm[3][1][0] = -tmp1 * xy - tmp3 - lx * glm[3].y - ly * glm[3].x;
    hlm[3][1][1] = -tmp1 * y2 - tmp0 - tmp2 - (ly + ly) * glm[3].y;
    hlm[3][0][2] = hlm[3][2][0] = -tmp1 * xz - lx * glm[3].z - lz * glm[3].x;
    hlm[3][1][2] = hlm[3][2][1] = -tmp1 * yz - ly * glm[3].z - lz * glm[3].y;
    hlm[3][2][2] = -tmp1 * z2 - tmp0 - (lz + lz) * glm[3].z;
    tmp0 = 1 * ylm[1];
    tmp1 = -1 * tmp0;
    hlm[1][0][0] = -tmp1 * x2 - tmp0 + tmp3 - (lx + lx) * glm[1].x;
    hlm[1][0][1] = hlm[1][1][0] = -tmp1 * xy + tmp2 - lx * glm[1].y - ly * glm[1].x;
    hlm[1][1][1] = -tmp1 * y2 - tmp0 - tmp3 - (ly + ly) * glm[1].y;
    hlm[1][0][2] = hlm[1][2][0] = -tmp1 * xz - lx * glm[1].z - lz * glm[1].x;
    hlm[1][1][2] = hlm[1][2][1] = -tmp1 * yz - ly * glm[1].z - lz * glm[1].y;
    hlm[1][2][2] = -tmp1 * z2 - tmp0 - (lz + lz) * glm[1].z;
    lx += x;
    ly += y;
    lz += z;
    tmp0 = qlm[7] * cm;
    tmp1 = qlm[7] * sm;
    tmp2 = 2.2360679774997898f * ylm[3];
    tmp3 = 2.2360679774997898f * ylm[1];
    glm[7].x = -lx * ylm[7] + tmp0;
    glm[7].y = -ly * ylm[7] - tmp1;
    glm[7].z = -lz * ylm[7] + tmp2;
    glm[5].x = -lx * ylm[5] + tmp1;
    glm[5].y = -ly * ylm[5] + tmp0;
    glm[5].z = -lz * ylm[5] + tmp3;
    tmp0 = 2 * ylm[7];
    tmp1 = 0 * tmp0;
    tmp2 = qlm[7] * cs;
    tmp3 = qlm[7] * ss;
    tmp4 = 2.2360679774997898f * qlm[4];
    tmp5 = tmp4 * sm;
    tmp4 *= cm;
    hlm[7][0][0] = -tmp1 * x2 - tmp0 + tmp2 - (lx + lx) * glm[7].x;
    hlm[7][0][1] = hlm[7][1][0] = -tmp1 * xy - tmp3 - lx * glm[7].y - ly * glm[7].x;
    hlm[7][1][1] = -tmp1 * y2 - tmp0 - tmp2 - (ly + ly) * glm[7].y;
    hlm[7][0][2] = hlm[7][2][0] = -tmp1 * xz + tmp4 - lx * glm[7].z - lz * glm[7].x;
    hlm[7][1][2] = hlm[7][2][1] = -tmp1 * yz - tmp5 - ly * glm[7].z - lz * glm[7].y;
    hlm[7][2][2] = -tmp1 * z2 - tmp0 - (lz + lz) * glm[7].z;
    tmp0 = 2 * ylm[5];
    tmp1 = 0 * tmp0;
    hlm[5][0][0] = -tmp1 * x2 - tmp0 + tmp3 - (lx + lx) * glm[5].x;
    hlm[5][0][1] = hlm[5][1][0] = -tmp1 * xy + tmp2 - lx * glm[5].y - ly * glm[5].x;
    hlm[5][1][1] = -tmp1 * y2 - tmp0 - tmp3 - (ly + ly) * glm[5].y;
    hlm[5][0][2] = hlm[5][2][0] = -tmp1 * xz + tmp5 - lx * glm[5].z - lz * glm[5].x;
    hlm[5][1][2] = hlm[5][2][1] = -tmp1 * yz + tmp4 - ly * glm[5].z - lz * glm[5].y;
    hlm[5][2][2] = -tmp1 * z2 - tmp0 - (lz + lz) * glm[5].z;
    cs = 2 * c1;
    ss = 2 * s1;
    cm = 2 * c0;
    sm = 2 * s0;
    s1 = s0;
    c1 = c0;
    c0 = x * c1 - y * s1;
    s0 = y * c1 + x * s1;
    lx = 2 * x;
    ly = 2 * y;
    lz = 2 * z;
    tmp0 = qlm[8] * cm;
    tmp1 = qlm[8] * sm;
    glm[8].x = -lx * ylm[8] + tmp0;
    glm[8].y = -ly * ylm[8] - tmp1;
    glm[8].z = -lz * ylm[8];
    glm[4].x = -lx * ylm[4] + tmp1;
    glm[4].y = -ly * ylm[4] + tmp0;
    glm[4].z = -lz * ylm[4];
    tmp0 = 2 * ylm[8];
    tmp1 = 0 * tmp0;
    tmp2 = qlm[8] * cs;
    tmp3 = qlm[8] * ss;
    hlm[8][0][0] = -tmp1 * x2 - tmp0 + tmp2 - (lx + lx) * glm[8].x;
    hlm[8][0][1] = hlm[8][1][0] = -tmp1 * xy - tmp3 - lx * glm[8].y - ly * glm[8].x;
    hlm[8][1][1] = -tmp1 * y2 - tmp0 - tmp2 - (ly + ly) * glm[8].y;
    hlm[8][0][2] = hlm[8][2][0] = -tmp1 * xz - lx * glm[8].z - lz * glm[8].x;
    hlm[8][1][2] = hlm[8][2][1] = -tmp1 * yz - ly * glm[8].z - lz * glm[8].y;
    hlm[8][2][2] = -tmp1 * z2 - tmp0 - (lz + lz) * glm[8].z;
    tmp0 = 2 * ylm[4];
    tmp1 = 0 * tmp0;
    hlm[4][0][0] = -tmp1 * x2 - tmp0 + tmp3 - (lx + lx) * glm[4].x;
    hlm[4][0][1] = hlm[4][1][0] = -tmp1 * xy + tmp2 - lx * glm[4].y - ly * glm[4].x;
    hlm[4][1][1] = -tmp1 * y2 - tmp0 - tmp3 - (ly + ly) * glm[4].y;
    hlm[4][0][2] = hlm[4][2][0] = -tmp1 * xz - lx * glm[4].z - lz * glm[4].x;
    hlm[4][1][2] = hlm[4][2][1] = -tmp1 * yz - ly * glm[4].z - lz * glm[4].y;
    hlm[4][2][2] = -tmp1 * z2 - tmp0 - (lz + lz) * glm[4].z;
}

void SHGradientL2(const float3& normDir, std::array<float, 9>& ylm, std::array<float3, 9>& glm)
{
    float x = normDir.x;
    float y = normDir.y;
    float z = normDir.z;
    float c0, c1, s0, s1, tmp, tmp0, tmp1, tmp2, tmp3;
    float z2 = z * z;
    std::array< float, 9 > qlm{};
    //zonal harmonics(m=0)
    ylm[0] = 0.2820947917738781f;
    qlm[0] = 0.2820947917738781f;
    ylm[2] = 0.4886025119029199f * z;
    qlm[3] = 0.4886025119029199f * z;
    ylm[6] = 0.9461746957575601f * z2 - 0.3153915652525200f;
    qlm[6] = 0.9461746957575601f * z2 - 0.3153915652525200f;
    c0 = x; s0 = y; c1 = 1; s1 = 0;
    //m = 001
    qlm[4] = -0.4886025119029200f;
    ylm[3] = qlm[4] * c0;
    ylm[1] = qlm[4] * s0;
    qlm[7] = -1.0925484305920792f * z;
    ylm[7] = qlm[7] * c0;
    ylm[5] = qlm[7] * s0;
    s1 = s0; c1 = c0; c0 = x * c1 - y * s1; s0 = y * c1 + x * s1;
    //m = 002
    qlm[8] = 0.5462742152960396f;
    ylm[8] = qlm[8] * c0;
    ylm[4] = qlm[8] * s0;
    //calculate gradient
    //const float3 zero = float3(0.f, 0.f, 0.f);
    glm[0] = float3(0.f, 0.f, 0.f);
    glm[2].x = -x * ylm[2];
    glm[2].y = -y * ylm[2];
    glm[2].z = -z * ylm[2] + 0.4886025119029199f;
    tmp0 = 2 * ylm[6]; tmp1 = 1.2909944487358056f * qlm[4];
    glm[6].x = -x * (tmp0 - tmp1);
    glm[6].y = -y * (tmp0 - tmp1);
    glm[6].z = -z * tmp0 + 2.5819888974716116f * qlm[3];
    c0 = x; s0 = y; c1 = 1; s1 = 0;
    //m = 001
    tmp = 1.2247448713915892f * qlm[2]; tmp0 = tmp * c0; tmp1 = tmp * s0;
    tmp = 1 * qlm[4]; tmp2 = tmp * c1; tmp3 = tmp * s1;
    tmp = 1 * ylm[3];
    glm[3].x = -x * (tmp - tmp0) + tmp2;
    glm[3].y = -y * (tmp - tmp0) - tmp3;
    glm[3].z = -z * tmp;
    tmp = 1 * ylm[1];
    glm[1].x = -x * (tmp - tmp1) + tmp3;
    glm[1].y = -y * (tmp - tmp1) + tmp2;
    glm[1].z = -z * tmp;
    tmp = 2.4494897427831783f * qlm[1]; glm[3].z += tmp * c0; glm[1].z += tmp * s0;
    tmp = 0.5270462766947299f * qlm[5]; tmp0 = tmp * c0; tmp1 = tmp * s0;
    tmp = 1 * qlm[7]; tmp2 = tmp * c1; tmp3 = tmp * s1;
    tmp = 2 * ylm[7];
    glm[7].x = -x * (tmp - tmp0) + tmp2;
    glm[7].y = -y * (tmp - tmp0) - tmp3;
    glm[7].z = -z * tmp;
    tmp = 2 * ylm[5];
    glm[5].x = -x * (tmp - tmp1) + tmp3;
    glm[5].y = -y * (tmp - tmp1) + tmp2;
    glm[5].z = -z * tmp;
    tmp = 2.2360679774997898f * qlm[4]; glm[7].z += tmp * c0; glm[5].z += tmp * s0;
    s1 = s0; c1 = c0; c0 = x * c1 - y * s1; s0 = y * c1 + x * s1;
    tmp = 2 * qlm[8]; tmp0 = tmp * c1; tmp1 = tmp * s1; tmp2 = 2 * ylm[8]; tmp3 = 2 * ylm[4];
    glm[8].x = -x * tmp2 + tmp0;
    glm[8].y = -y * tmp2 - tmp1;
    glm[8].z = -z * tmp2;
    glm[4].x = -x * tmp3 + tmp1;
    glm[4].y = -y * tmp3 + tmp0;
    glm[4].z = -z * tmp3;
}


//void decomposeSH(std::vector<float4>& out, const Falcor::ref<EnvMap>& envMap)
//{
//    if (shOrder == -1)
//    {
//        logError( "call initSHTable before decompositionSHEnvMap!");
//        return;
//    }
//
//    int num_basis = (shOrder + 1) * (shOrder + 1);
//    out.resize(num_basis);
//    const Falcor::ref<Texture> envMapTex = envMap->getEnvMap();
//    int width = envMapTex->getWidth();
//    int height = envMapTex->getHeight();
//
//    //read texture data into an array
//    float4* data = new float4[width*height];
//    envMapTex->getSubresourceBlob(0, &data[0], sizeof(float4) * width * height);
//
//   // float4* tranposedData = TranposeData(data, width, height);
//
//    for (int l = 0; l < num_basis; ++l)
//    {
//        double r = 0.0;
//        double g = 0.0;
//        double b = 0.0;
//        double a = 0.0;
//
//        for (int y = 0; y < height; ++y)
//        {
//            for (int x = 0; x < width; ++x)
//            {
//                //read value from env map
//                float4 envMapValue = data[TWO_D_TO_ONE_D(x,y,width)]; 
//
//                // Lookup SH basis value at (l, x, y)
//                float yd = SHBasisTable[THREE_D_TO_ONE_D(l, x, y, width, height)]; 
//
//               //  precomputed pixel solid angle weight sin(theta)*dtheta*dphi
//                float delta = dOmega[TWO_D_TO_ONE_D(x, y, width)]; 
//
//                double weight = (double)yd * (double)delta;
//
//                r += (double)envMapValue[0] * weight;
//                g += (double)envMapValue[1] * weight;
//                b += (double)envMapValue[2] * weight;
//                a += (double)envMapValue[3] * weight;
//            }
//        }
//        float4 tmp(r, g, b, a);
//        out[l] = tmp;
//    }
//    delete [] data;
//}
//
//void reconstructSH(std::vector<float4>& sh_coeff, const Falcor::ref<EnvMap>& envMap, Falcor::ref<Device> pDevice)
//{
//    int num_basis = (shOrder + 1) * (shOrder + 1);
//    const Falcor::ref<Texture> envMapTex = envMap->getEnvMap();
//    int width = envMapTex->getWidth();
//    int height = envMapTex->getHeight();
//
//  // float3* reconstructedData = new float3[width * height];
//
//    //for (int y = 0; y < height; ++y)
//    //{
//    //    for (int x = 0; x < width; ++x)
//    //    {
//    //        float4 tmp = float4 (0,0,0,0);
//    //        for (int i = 0; i < num_basis; ++i)
//    //        {
//    //            tmp += sh_coeff[i] * envSHTable[THREE_D_TO_ONE_D(i, x, y, width, height)]; 
//    //        }
//    //        float3 color = float3(tmp[0], tmp[1], tmp[2]);
//    //        reconstructedData[TWO_D_TO_ONE_D(x, y, width)] = color;
//    //    }
//    //}
//
//   float* reconstructedData = new float[width * height*3];
//
//    for (int y = 0; y < height; ++y)
//    {
//        for (int x = 0; x < width; ++x)
//        {
//            float4 tmp = float4 (0,0,0,0);
//            for (int i = 0; i < num_basis; ++i)
//            {
//                tmp += sh_coeff[i] * SHBasisTable[THREE_D_TO_ONE_D(i, x, y, width, height)]; 
//            }
//            float3 color = float3(tmp[0], tmp[1], tmp[2]);
//            reconstructedData[TWO_D_TO_ONE_D(x, y, width)*3 + 0] = color[0];
//            reconstructedData[TWO_D_TO_ONE_D(x, y, width)*3 + 1] = color[1];
//            reconstructedData[TWO_D_TO_ONE_D(x, y, width)*3 + 2] = color[2];
//        }
//    }

    //Falcor::ref<Texture> outTex =
    //    pDevice->createTexture2D(width, height, ResourceFormat::RGB32Float, 1, 1, reconstructedData, ResourceBindFlags::ShaderResource);

    //outTex->captureToFile(0, 0, "reconstructed.pfm", Bitmap::FileFormat::PfmFile, Bitmap::ExportFlags::None, false);
    //cout << "First 3" << endl;
    //for (int i = 0; i < 9 ; i+=3)
    //{
    //    cout << reconstructedData[i] << " " << reconstructedData[i + 1] << " " << reconstructedData[i + 2] << endl;
    //}
    //cout << "end first 3" << endl;

//    stbi_write_hdr("reconstructed_env.hdr", width, height, 3, reconstructedData);
//
//    delete[] reconstructedData;
//}

float4* TranposeData(float4* data, int width, int height)
{
    // tranpose data because input hdr file is column major
    float4* tranposedData = new float4[width * height];

    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            int srcIdx = (x * height + y); // column-major
            int dstIdx = (y * width + x);  // row-major

            tranposedData[dstIdx + 0] = data[srcIdx + 0];
        }
    }
    return tranposedData;
}

void computeProbesPos(ProbeGrid& grid)
{
    const int3 res = grid.resolution;
    const float3 origin = grid.origin;
    const float3 spacing = grid.spacing;
    int numProbes = grid.resolution.x * grid.resolution.y * grid.resolution.z;

    int width = res.x;
    int height = res.y;
    int depth = res.z;

    int num_basis = grid.numBasis;
    grid.probesSHCoeffs.resize(numProbes * num_basis);
    
    grid.probesPos.reserve(numProbes);
    for (int probeZ = 0; probeZ < depth; ++probeZ)
    {
        for (int probeY = 0; probeY < height; ++probeY)
        {
            for (int probeX = 0; probeX < width; ++probeX)
            {
                //int index = probeX + probeY * width + probeZ * width * height;
                float3 probePos = {origin.x + probeX * spacing.x, origin.y + probeY * spacing.y, origin.z + probeZ * spacing.z};
                grid.probesPos.push_back(probePos);
            }
        }
    }
}

void saveProbeGridToFile(const ProbeGrid& grid, const std::string& path)
{
    std::ofstream file(path, std::ios::binary);
    if (!file)
        return;

    // Set high precision for floating point values
    file << std::fixed << std::setprecision(6);
    int numBasis = grid.numBasis;
    file << "# Probe Grid Info\n";
    file << "Resolution: " << grid.resolution.x << " " << grid.resolution.y << " " << grid.resolution.z << "\n";
    file << "Origin: " << grid.origin.x << " " << grid.origin.y << " " << grid.origin.z << "\n";
    file << "Spacing: " << grid.spacing.x << " " << grid.spacing.y << " " << grid.spacing.z << "\n";
    file << "NumBasis: " << numBasis << "\n";

    file << "# Probes (SH Coefficients per Probe)\n";

    int probeCount = 0;
    int numProbes = grid.resolution.x * grid.resolution.y * grid.resolution.z;

    for (int p = 0; p < numProbes; ++p)
    {
        file << "Probe " << probeCount++ << ":\n";
        file << "Coeffs:\n";
        for (int b = 0; b < numBasis; ++b)
        {
            const float4& coeff = grid.probesSHCoeffs[p * numBasis + b];
            file << "  " << coeff.x << " " << coeff.y << " " << coeff.z << " " << coeff.w << "\n";
        }
        file << "\n";
    }

    file.close();
}

void saveProbeGridToFileWithGradAndHessian(const ProbeGrid& grid, const std::string& path)
{
    std::ofstream file(path, std::ios::binary);
    if (!file)
        return;

    // Set high precision for floating point values
    file << std::fixed << std::setprecision(6);
    int numBasis = grid.numBasis;
    file << "# Probe Grid Info\n";
    file << "Resolution: " << grid.resolution.x << " " << grid.resolution.y << " " << grid.resolution.z << "\n";
    file << "Origin: " << grid.origin.x << " " << grid.origin.y << " " << grid.origin.z << "\n";
    file << "Spacing: " << grid.spacing.x << " " << grid.spacing.y << " " << grid.spacing.z << "\n";
    file << "NumBasis: " << numBasis << "\n";

    file << "# Probes (SH Coefficients per Probe)\n";

    int probeCount = 0;
    int numProbes = grid.resolution.x * grid.resolution.y * grid.resolution.z;

    for (int probeIdx = 0; probeIdx < numProbes; ++probeIdx)
    {
        file << "Probe " << probeCount++ << ":\n";
        file << "Coeffs:\n";
        for (int basisIdx = 0; basisIdx < numBasis; ++basisIdx)
        {
            const float4& coeff = grid.probesSHCoeffs[probeIdx * numBasis + basisIdx];
            file << "  " << coeff.x << " " << coeff.y << " " << coeff.z << " " << coeff.w << "\n";
        }
        file << "\n";
        file << "*** Gradients: (each line corresponds to r, g, b channel respectively) ***\n";
        for (int basisIdx = 0; basisIdx < numBasis; ++basisIdx)
        {
            file << "*** of coeff " << basisIdx << " ***\n";
            const GradSHCoeff& coeff = grid.probesSHCoeffsGradients[probeIdx * numBasis + basisIdx];
            file << coeff.r.x << " " << coeff.r.y << " " << coeff.r.z << "\n";
            file << coeff.g.x << " " << coeff.g.y << " " << coeff.g.z << "\n";
            file << coeff.b.x << " " << coeff.b.y << " " << coeff.b.z << "\n";
            file << "\n";
        }
    
        file << "*** Hessian: (each 3 lines corresponds to r, g, b channel respectively) ***\n";
        for (int basisIdx = 0; basisIdx < numBasis; ++basisIdx)
        {
            file << "*** of coeff " << basisIdx << " ***\n";
            const HessianSHCoeff& hessian = grid.probesSHCoeffsHessians[probeIdx * numBasis + basisIdx];
            for (int row = 0; row < 3; ++row)
            {
                file << hessian.r[row][0] << " " << hessian.r[row][1] << " " << hessian.r[row][2] << "\n";
            }
            file << "\n";
            for (int row = 0; row < 3; ++row)
            {
                file << hessian.g[row][0] << " " << hessian.g[row][1] << " " << hessian.g[row][2] << "\n";
            }
            file << "\n";
            for (int row = 0; row < 3; ++row)
            {
                file << hessian.b[row][0] << " " << hessian.b[row][1] << " " << hessian.b[row][2] << "\n";
            }
            file << "\n";
        }
    }

    file.close();
}

bool loadProbeGridFromFileWithGradAndHessian(ProbeGrid& grid, const std::string& path)
{
    std::ifstream file(path);
    if (!file)
        return false;

    std::string line;
    int numBasis = 0;

    std::getline(file, line); // Skip "# Probe Grid Info"

    // Resolution
    std::getline(file, line);
    std::istringstream resStream(line.substr(line.find(":") + 1));
    resStream >> grid.resolution.x >> grid.resolution.y >> grid.resolution.z;

    // Origin
    std::getline(file, line);
    std::istringstream orgStream(line.substr(line.find(":") + 1));
    orgStream >> grid.origin.x >> grid.origin.y >> grid.origin.z;

    // Spacing
    std::getline(file, line);
    std::istringstream spcStream(line.substr(line.find(":") + 1));
    spcStream >> grid.spacing.x >> grid.spacing.y >> grid.spacing.z;

    // NumBasis
    std::getline(file, line);
    std::istringstream basisStream(line.substr(line.find(":") + 1));
    basisStream >> numBasis;
    grid.numBasis = numBasis;
    
    std::getline(file, line); // Skip "# Probes (SH Coefficients per Probe)"

    grid.probesSHCoeffs.clear();
    int numProbes = grid.resolution.x * grid.resolution.y * grid.resolution.z;
    grid.probesSHCoeffs.resize(numProbes * numBasis);

    grid.probesSHCoeffsGradients.clear();
    grid.probesSHCoeffsGradients.resize(numProbes * numBasis);

    grid.probesSHCoeffsHessians.clear();
    grid.probesSHCoeffsHessians.resize(numProbes * numBasis);

    for (int probeIdx = 0; probeIdx < numProbes; ++probeIdx)
    {
        // Read SH coefficients
        std::getline(file, line); // skip "Probe X"
        std::getline(file, line); // skip "Coeffs:"

        for (int basisIdx = 0; basisIdx < numBasis; ++basisIdx)
        {
            std::getline(file, line); // get coeff line
            std::istringstream coeffStream(line);
            float4 coeff;
            coeffStream >> coeff.x >> coeff.y >> coeff.z >> coeff.w;
            grid.probesSHCoeffs[probeIdx * numBasis + basisIdx] = coeff;
        }

         // Read gradients
        std::getline(file, line); // empty line
        std::getline(file, line); // skip "*** Gradients: (each line corresponds to r, g, b channel respectively) ***"
        
        for (int basisIdx = 0; basisIdx < numBasis; ++basisIdx)
        {
            std::getline(file, line); // skip " *** of coeff X ***"
            float3 gradR, gradG, gradB;
            for (int channel = 0; channel < 3; ++channel)
            {
                std::getline(file, line); // get gradient lines
                std::istringstream gradStream(line);
                if (channel == 0)
                    gradStream >> gradR.x >> gradR.y >> gradR.z;
                else if (channel == 1)
                    gradStream >> gradG.x >> gradG.y >> gradG.z;
                else
                    gradStream >> gradB.x >> gradB.y >> gradB.z;
            }
            grid.probesSHCoeffsGradients[probeIdx * numBasis + basisIdx] = {gradR, gradG, gradB};
            std::getline(file, line); // skip empty line
        }

        //Read Hessians
        std::getline(file, line); // skip "*** Hessian: (each 3 lines corresponds to r, g, b channel respectively) ***"
        for (int basisIdx = 0; basisIdx < numBasis; ++basisIdx)
        {
            std::getline(file, line); // skip " *** of coeff X ***"
            float3x3 hessR, hessG, hessB;
            for (int channel = 0; channel < 3; ++channel)
            {
                for (int row = 0; row < 3; ++row)
                {
                    std::getline(file, line); // get hessian lines
                    std::istringstream hessStream(line);
                    if (channel == 0)
                        hessStream >> hessR[row][0] >> hessR[row][1] >> hessR[row][2];
                    else if (channel == 1)
                        hessStream >> hessG[row][0] >> hessG[row][1] >> hessG[row][2];
                    else
                        hessStream >> hessB[row][0] >> hessB[row][1] >> hessB[row][2];
                }
                std::getline(file, line); // skip empty line
            }
            grid.probesSHCoeffsHessians[probeIdx * numBasis + basisIdx] = {hessR, hessG, hessB};
        }
    }

    computeProbesPos(grid);
    file.close();
    return grid.probesSHCoeffs.size() == numProbes;
}

bool loadProbeGridFromFile(ProbeGrid& grid, const std::string& path)
{
    std::ifstream file(path);
    if (!file)
        return false;

    std::string line;
    int numBasis = 0;

    std::getline(file, line); // Skip "# Probe Grid Info"

    // Resolution
    std::getline(file, line);
    std::istringstream resStream(line.substr(line.find(":") + 1));
    resStream >> grid.resolution.x >> grid.resolution.y >> grid.resolution.z;

    // Origin
    std::getline(file, line);
    std::istringstream orgStream(line.substr(line.find(":") + 1));
    orgStream >> grid.origin.x >> grid.origin.y >> grid.origin.z;

    // Spacing
    std::getline(file, line);
    std::istringstream spcStream(line.substr(line.find(":") + 1));
    spcStream >> grid.spacing.x >> grid.spacing.y >> grid.spacing.z;

    // NumBasis
    std::getline(file, line);
    std::istringstream basisStream(line.substr(line.find(":") + 1));
    basisStream >> numBasis;
    grid.numBasis = numBasis;
    // Skip "# Probes (SH Coefficients per Probe)"
    std::getline(file, line);

    grid.probesSHCoeffs.clear();
    int numProbes = grid.resolution.x * grid.resolution.y * grid.resolution.z;
    grid.probesSHCoeffs.resize(numProbes*numBasis);

    for (int probeIndex = 0; probeIndex < numProbes; ++probeIndex)
    {
        std::getline(file, line); // "Probe X"

        for (int i = 0; i < numBasis; ++i)
        {
            if (!std::getline(file, line))
                return false;
            std::istringstream coeffStream(line);
            float4 coeff;
            coeffStream >> coeff.x >> coeff.y >> coeff.z >> coeff.w;
            grid.probesSHCoeffs[probeIndex * numBasis + i] = coeff;
        }

        std::getline(file, line); // empty line
    }

    computeProbesPos(grid);
    return grid.probesSHCoeffs.size() == numProbes;
}

void generateUniformSphereDirSamples(int sampleCount, std::vector<ProbeDirSample>& out)
{
    out.clear();
    out.reserve(sampleCount);

    //std::mt19937 rng(12345); // fixed seed for reproducibility
    //std::mt19937 rng(123); // fixed seed for reproducibility
    //std::uniform_real_distribution<float> dist(0.0f, 1.0f);

    //float dOmega = 4.0f * float(M_PI) / float(sampleCount); // Uniform solid angle per sample

    //for (int i = 0; i < sampleCount; ++i)
    //{
    //     Generate in standard polar coordinates z -up, y -right, x -forward
    //    float z = 1.0f - 2.0f * dist(rng);          // up 
    //    float phi = 2.0f * float(M_PI) * dist(rng); // azimuth
    //    float h = sqrtf(1.0f - z * z);              // radius in xy-plane

    //    
    //    float x = h * cosf(phi); // forward
    //    float y = h * sinf(phi); // right

    //    float3 dir = {x, y, z};  // x, y, z
    //    dir = math::normalize(dir);
    //    out.push_back({dir, dOmega});
    //}

    float phi = (sqrtf(5.0f) - 1.0f) / 2.0f; // Golden Ratio
    float ga = phi * 2.0f * float(M_PI);     // Golden Angle

    float dOmega = 4.0f * float(M_PI) / float(sampleCount);

    for (int i = 0; i < sampleCount; ++i)
    {
        // 1. Z goes from 1 to -1 evenly
        // Adding 0.5 to i aligns samples better at poles
        float z = (1.0f - (2.0f * float(i) + 1.0f) / float(sampleCount));

        // 2. Radius at this Z
        float radius = sqrtf(1.0f - z * z);

        // 3. Golden Angle increment
        float theta = ga * float(i);

        float x = radius * cosf(theta);
        float y = radius * sinf(theta);

        float3 dir = { x, y, z };
        dir = math::normalize(dir); // robust normalization
        out.push_back({ dir, dOmega });
    }
}

float3 gradientOmega(float3 s, float3 x, float3 n, float N)
{
    // q = s - x
    float3 q = s - x;

    // r = ||q||
    float r = length(q);

    // cosξ = -(n · q) / r
    float cosXi = -(dot(n, q)) / r;

    // Numerical guard to avoid division by zero or near-zero values
    const float eps = 1e-6f;
    if (fabs(cosXi) < eps)
        cosXi = (cosXi < 0 ? -eps : eps);

    // Use uniform sampling with N samples
    float factor = 4.0f * M_PI / float(N);

    // Element-wise computation:
    // ∂xΩ_i = factor * (n_x * r + 3 * q_x * cosξ) / (r² * cosξ)
    // ∂yΩ_i = factor * (n_y * r + 3 * q_y * cosξ) / (r² * cosξ)
    // ∂zΩ_i = factor * (n_z * r + 3 * q_z * cosξ) / (r² * cosξ)
    float3 grad = factor * ((n * r + 3.0f * q * cosXi) / (r * r * cosXi));

    return grad;
}



float3x3 hessianOmega(const float3& s, const float3& x, const float3& n, int N)
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

    float factor = 4.0f * M_PI / float(N);

    // Pure second derivatives
    float d2Omega_xx = factor * (6.0f * n.x * qx * r - 3.0f * cosXi * (r2 - 5.0f * qx * qx)) / (r4 * cosXi);
    float d2Omega_yy = factor * (6.0f * n.y * qy * r - 3.0f * cosXi * (r2 - 5.0f * qy * qy)) / (r4 * cosXi);
    float d2Omega_zz = factor * (6.0f * n.z * qz * r - 3.0f * cosXi * (r2 - 5.0f * qz * qz)) / (r4 * cosXi);

    // Mixed derivatives (symmetric)
    float d2Omega_xy = factor * (((3.0f * n.x * qy + 3.0f * qx * n.y) / (r3 * cosXi)) + (15.0f * qx * qy / r4));
    float d2Omega_xz = factor * (((3.0f * n.x * qz + 3.0f * qx * n.z) / (r3 * cosXi)) + (15.0f * qx * qz / r4));
    float d2Omega_yz = factor * (((3.0f * n.y * qz + 3.0f * qy * n.z) / (r3 * cosXi)) + (15.0f * qy * qz / r4));

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

// Notes:
//   - gradOmega() computes ∇Ω_i (geometry term derivative).
//   - gradY_lm() provides ∇Y_l^m(ω_i) in direction space.
//   - ∂ω/∂x = -(I - ωω^T) / r maps direction derivative to spatial domain.
//   - Numerical guards in gradOmega() handle small r or cosξ cases.
//   - Typically used for SH grid refinement or Hessian computation.
//-------------------------------------------------------------------------------
void calculateGradAndHessianSHCoeffLM(
    const float3& x,
    const std::vector<ProbeSampleData>& samplingData,
    const std::vector<ProbeDirSample>& samplingDir,
    const int& basisIdx,
    GradSHCoeff& outGrad,
    HessianSHCoeff& outHessian
)
{
    outGrad = {float3(0), float3(0), float3(0)};
    outHessian = {float3x3::zeros(), float3x3::zeros(), float3x3::zeros()};
    int samplingSize = samplingData.size();
    for (int sampleIdx = 0; sampleIdx < samplingSize; ++sampleIdx)
    {
        if (samplingData[sampleIdx].hitT < 0.0f) // ray miss
            continue;
        float3 s = float3(samplingData[sampleIdx].s.x, samplingData[sampleIdx].s.y, samplingData[sampleIdx].s.z);
        float3 n = float3(samplingData[sampleIdx].n.x, samplingData[sampleIdx].n.y, samplingData[sampleIdx].n.z);
        float3 L = float3(samplingData[sampleIdx].Li.x, samplingData[sampleIdx].Li.y, samplingData[sampleIdx].Li.z);

        // Solid angle (uniform per patch)
        float Omega_i = (4.0f * M_PI) / (float)samplingSize;

        // Gradient ∂_x Ω_i
        float3 gradOmega = gradientOmega(s, x, n, samplingSize);

        // SH basis and derivative in direction space
        int numBasis = 9;
        float Ylm = SHBasisTable[numBasis * sampleIdx + basisIdx];
        float3 gradYlm = SHGradientTable[numBasis * sampleIdx + basisIdx];

        // q = s - x
        float3 q = s - x;

        // r = ||q||
        float rInv = 1.0f / length(q);

        // Contribution of this sample to ∇f_l^m
        // ∇_x Ω_i Y_l^m (ω_i )-Ω_i/r_i  ∇_(ω_i ) Y_l^m (ω_i)
        float3 contrib = gradOmega * Ylm - Omega_i * rInv * gradYlm;
        outGrad.r += (L.r * contrib);
        outGrad.g += (L.g * contrib);
        outGrad.b += (L.b * contrib);

        float3x3 H_Omega = hessianOmega(s, x, n, samplingSize);
        float3x3 hessYlm = SHHessianTable[numBasis * sampleIdx + basisIdx];

        // Pre-calculate shared values for the inner loop
        float3 w = samplingDir[sampleIdx].dir;
        float rInvSq = rInv * rInv;

        //  Accumulate Hessian of f_l^m
        // ∂_(x_j x_k)
        //    f_l ^m =∑_(i = 1) ^ N▒L(ω_i)[(∂_(x_j x_k) Ω_i)Y_l ^ m - 1 / r_i((∂_(x_j) Ω_i)(∂_(ω_i, k) Y_l ^ m) + (∂_(x_k) Ω_i)(∂_(ω_i, j)
        //    Y_l ^ m)) +
        //                                                Ω_i / (r_i ^ 2)((ω_i)_j(∂_(ω_i, k) Y_l ^ m) + ∂_(ω_i, j) ∂_(ω_i, k) Y_l ^ m)]
        for (int j = 0; j < 3; ++j)
        {
            for (int k = 0; k < 3; ++k)
            {
                // Term 1: (∂²Ω / ∂xj∂xk) * Y
                float term1 = H_Omega[j][k] * Ylm;

                // Term 2: -1/r * [ (∂Ω/∂xj)(∂Y/∂wk) + (∂Ω/∂xk)(∂Y/∂wj) ]
                // Note: We access float3 components using array indexing [j] and [k]
                float term2 = -rInv * (gradOmega[j] * gradYlm[k] + gradOmega[k] * gradYlm[j]);

                // Term 3: Ω/r² * [ (ω)_j * (∂Y/∂wk) + (∂²Y / ∂wj∂wk) ]
                float term3 = Omega_i * rInvSq * (w[j] * gradYlm[k] + hessYlm[j][k]);

                float contrib = term1 + term2 + term3;
                outHessian.r[j][k] += L.r * contrib;
                outHessian.g[j][k] += L.g * contrib;
                outHessian.b[j][k] += L.b * contrib;
            }
        }
    }
}

void calculateGradSHCoeffLM(
    const float3& x,
    const std::vector<ProbeSampleData>& samplingData,
    const std::vector<ProbeDirSample>& samplingDir,
    const int& basisIdx,
    GradSHCoeff& outGrad
)
{
    outGrad = { float3(0), float3(0), float3(0) };
    int samplingSize = samplingData.size();
    for (int sampleIdx = 0; sampleIdx < samplingSize; ++sampleIdx)
    {
        if (samplingData[sampleIdx].hitT < 0.0f) // ray miss
            continue;
        float3 s = float3(samplingData[sampleIdx].s.x, samplingData[sampleIdx].s.y, samplingData[sampleIdx].s.z);
        float3 n = float3(samplingData[sampleIdx].n.x, samplingData[sampleIdx].n.y, samplingData[sampleIdx].n.z);
        float3 L = float3(samplingData[sampleIdx].Li.x, samplingData[sampleIdx].Li.y, samplingData[sampleIdx].Li.z);

        // Solid angle (uniform per patch)
        float Omega_i = (4.0f * M_PI) / (float)samplingSize;

        // Gradient ∂_x Ω_i
        float3 gradOmega = gradientOmega(s, x, n, samplingSize);

        // SH basis and derivative in direction space
        int numBasis = 9;
        float Ylm = SHBasisTable[numBasis * sampleIdx + basisIdx];
        float3 gradYlm = SHGradientTable[numBasis * sampleIdx + basisIdx];

        // q = s - x
        float3 q = s - x;

        // r = ||q||
        float rInv = 1.0f / length(q);

        // Contribution of this sample to ∇f_l^m
        // ∇_x Ω_i Y_l^m (ω_i )-Ω_i/r_i  ∇_(ω_i ) Y_l^m (ω_i)
        float3 contrib = gradOmega * Ylm - Omega_i * rInv * gradYlm;
        outGrad.r += (L.r * contrib);
        outGrad.g += (L.g * contrib);
        outGrad.b += (L.b * contrib);
    }
}

HessianSHCoeff hessianSHCoeffLM(const float3& x, const std::vector<ProbeSampleData>& samplingData, const int& basisIdx)
{
    HessianSHCoeff H = {float3x3::zeros(), float3x3::zeros(), float3x3::zeros()};

    int numBasis = 9;
    int samplingSize = samplingData.size();
    for (int sampleIdx = 0; sampleIdx < samplingSize; ++sampleIdx)
    {
        if (samplingData[sampleIdx].hitT < 0.0f) // ray miss
            continue;
        float3 s = float3(samplingData[sampleIdx].s.x, samplingData[sampleIdx].s.y, samplingData[sampleIdx].s.z);
        float3 n = float3(samplingData[sampleIdx].n.x, samplingData[sampleIdx].n.y, samplingData[sampleIdx].n.z);
        float3 L = float3(samplingData[sampleIdx].Li.x, samplingData[sampleIdx].Li.y, samplingData[sampleIdx].Li.z);

        float3 q = s - x;
        float r = length(q);
        float Omega_i = (4.0f * M_PI) / samplingSize;

        // Gradients and Hessian of Omega_i
        float3 gOmega = gradientOmega(s, x, n, samplingSize);
        float3x3 H_Omega = hessianOmega(s, x, n, samplingSize);

         // SH basis and derivatives
        float Ylm = SHBasisTable[numBasis * sampleIdx + basisIdx];
        float3 gradYlm = SHGradientTable[numBasis * sampleIdx + basisIdx];
        float3x3 hessYlm = SHHessianTable[numBasis * sampleIdx + basisIdx];

        // Accumulate Hessian of f_l^m
        //∂_(x_j x_k ) f_l^m=∑_(i=1)^N▒L(ω_i)((∂_(x_j x_k ) Ω_i)Y_l^m (ω_i)+(∂_(x_j ) Ω_i)(∂_(x_k ) Y_l^m (ω_i))+(∂_(x_k ) Ω_i)(∂_(x_j ) Y_l^m (ω_i))+Ω_i 〖(∂〗_(x_j x_k ) Y_l^m (ω_i)))
        for (int j = 0; j < 3; ++j)
        {
            for (int k = 0; k < 3; ++k)
            {
                float contrib = H_Omega[j][k] * Ylm + gOmega[j] * gradYlm[k] + gOmega[k] * gradYlm[j] + Omega_i * hessYlm[j][k];
                H.r[j][k] += L.x * contrib;
                H.g[j][k] += L.y * contrib;
                H.b[j][k] += L.z * contrib;
               // HessianCoeffLM[j][k] += L * (H_Omega[j][k] * Ylm + gOmega[j] * gradYlm[k] + gOmega[k] * gradYlm[j] + Omega_i * hessYlm[j][k]);
            }
        }
    }

    return H;
}

std::vector<float3> generateVerificationPositions(float y, float extent, uint32_t resolution, float h)
{
    std::vector<float3> positions;
    positions.reserve(resolution * resolution * 3);

    double step = (2.0 * extent) / (double)(resolution - 1);

    //for (uint32_t i = 0; i < resolution; ++i)
    //{
    //    double z = -extent + step * (double)i;

    //    for (uint32_t j = 0; j < resolution; ++j)
    //    {
    //        double x = -extent + step * (double)j;

    //        positions.push_back(float3((float)x, y, (float)z));       // center
    //        positions.push_back(float3((float)(x + h), y, (float)z)); // +h
    //        positions.push_back(float3((float)(x - h), y, (float)z)); // -h
    //    }
    //}

    float z = 0;
    for (uint32_t j = 0; j < resolution; ++j)
    {
        double x = -extent + step * (double)j;

        positions.push_back(float3((float)x, y, z));       // center
        positions.push_back(float3((float)(x + h), y, z)); // +h
        positions.push_back(float3((float)(x - h), y, z)); // -h
    }

    return positions;
}

std::vector<float3> generateVerificationPositionsMixed(float y, float extent, uint32_t resolution, float h)
{
    // Order per sample: [center, x+h z+h, x+h z-h, x-h z+h, x-h z-h]
    // Total 5 points per sample location
    std::vector<float3> positions;
    positions.reserve(resolution * 5);

    double step = (2.0 * extent) / (double)(resolution - 1);

    // Base Z is 0.0f
    float z_base = -0.25f;
    //float z_base = 0.0f;

    for (uint32_t j = 0; j < resolution; ++j)
    {
        double x_d = -extent + step * (double)j;
        float x = (float)x_d;

        // 1. Center
        positions.push_back(float3(x, y, z_base));

        // 2. Corner: x+h, z+h
        positions.push_back(float3(x + h, y, z_base + h));

        // 3. Corner: x+h, z-h
        positions.push_back(float3(x + h, y, z_base - h));

        // 4. Corner: x-h, z+h
        positions.push_back(float3(x - h, y, z_base + h));

        // 5. Corner: x-h, z-h
        positions.push_back(float3(x - h, y, z_base - h));
    }

    return positions;
}

float calculateChannelRSHCoeffLM(int basisIdx, const std::vector<ProbeSampleData>& probeSamplingResults)
{
    if (shOrder == -1 || SHBasisTable == nullptr)
    {
        logError("call initSHTable before calculate SH coeffs!");
        return -1;
    }

    int num_basis = 9;
    int numSamplePerProbe = probeSamplingResults.size();
    float weight = 4.0f * float(M_PI) / float(numSamplePerProbe); // because uniform sampling

    float coeffR = 0.0f;
    // For each direction sample
    for (int sampleIdx = 0; sampleIdx < numSamplePerProbe; ++sampleIdx)
    {
        const float4& sample = probeSamplingResults[sampleIdx].Li;
        float shY = SHBasisTable[num_basis * sampleIdx + basisIdx]; // SH basis value for this direction

        coeffR += sample.r * shY * weight;
    }

    return coeffR;
}

void calculateChannelRGradAndHessianSHCoeffLM(
    const float3& x,
    const std::vector<ProbeSampleData>& samplingData,
    const std::vector<ProbeDirSample>& samplingDir,
    const int& basisIdx,
    float3& outGrad,
    float3x3& outHessian
)
{
    outGrad = float3(0.0f, 0.0f, 0.0f);
    outHessian = float3x3::zeros();
    int samplingSize = samplingData.size();

    for (int sampleIdx = 0; sampleIdx < samplingSize; ++sampleIdx)
    {
        if (samplingData[sampleIdx].hitT < 0.0f) // ray miss
            continue;

        float3 s = float3(samplingData[sampleIdx].s.x, samplingData[sampleIdx].s.y, samplingData[sampleIdx].s.z);
        float3 n = float3(samplingData[sampleIdx].n.x, samplingData[sampleIdx].n.y, samplingData[sampleIdx].n.z);
        float3 L = float3(samplingData[sampleIdx].Li.x, samplingData[sampleIdx].Li.y, samplingData[sampleIdx].Li.z);

        // Solid angle (uniform per patch)
        float Omega_i = (4.0f * M_PI) / (float)samplingSize;

        // Gradient ∂_x Ω_i
        float3 gradOmega = gradientOmega(s, x, n, samplingSize);

        // SH basis and derivative in direction space
        int numBasis = 9;
        float Ylm = SHBasisTable[numBasis * sampleIdx + basisIdx];
        float3 gradYlm = SHGradientTable[numBasis * sampleIdx + basisIdx];

        // q = s - x
        float3 q = s - x;

        // r = ||q||
        float rInv = 1.0f / length(q);

        // Contribution of this sample to ∇f_l^m
        // ∇_x Ω_i Y_l^m (ω_i )-Ω_i/r_i  ∇_(ω_i ) Y_l^m (ω_i)
        float3 contrib = gradOmega * Ylm - Omega_i * rInv * gradYlm;
        outGrad += (L.r * contrib);

        float3x3 H_Omega = hessianOmega(s, x, n, samplingSize);
        float3x3 hessYlm = SHHessianTable[numBasis * sampleIdx + basisIdx];

        // Pre-calculate shared values for the inner loop
        float3 w = samplingDir[sampleIdx].dir;
        float rInvSq = rInv * rInv;

        // Accumulate Hessian of f_l^m
        for (int j = 0; j < 3; ++j)
        {
            for (int k = 0; k < 3; ++k)
            {
                // Term 1: (∂²Ω / ∂xj∂xk) * Y
                float term1 = H_Omega[j][k] * Ylm;

                // Term 2: -1/r * [ (∂Ω/∂xj)(∂Y/∂wk) + (∂Ω/∂xk)(∂Y/∂wj) ]
                float term2 = -rInv * (gradOmega[j] * gradYlm[k] + gradOmega[k] * gradYlm[j]);

                // Term 3: Ω/r² * [ (∂²Y / ∂wj∂wk) + (ω)_j * (∂Y/∂wk) + (ω)_k * (∂Y/∂wj) ]

                float term3 = Omega_i * rInvSq * hessYlm[j][k];

                float contrib = term1 + term2 + term3;
                outHessian[j][k] += L.r * contrib;
            }
        }
    }
}

float3 testComputeSHGrad()
{
     float3 testS = float3(1.0f, 0.0f, 0.0f);
     float3 testX = float3(0.0f, 0.0f, 0.0f);
     float3 testNormal = float3(-1.0f, 0.0f, 0.0f);
     int N = 4096;
     float3 gradOmega = gradientOmega(testS, testX, testNormal, N);

       // SH basis and derivative in direction space
     int numBasis = 9;
     /*float Ylm = SHBasisTable[numBasis * sampleIdx + basisIdx];
     float3 gradYlm = SHGradientTable[numBasis * sampleIdx + basisIdx];*/
     float Ylm = 0.282094792f;
     float L = 2.0f;
     float3 gradYlm = float3(0.0f, 0.0f, 0.0f); // for l=0, m=0, gradient is zero
    float Omega_i = (4.0f * M_PI) / (float)N;
     float3 contrib = gradOmega * Ylm + Omega_i * gradYlm;
    return L*contrib;
}

// Add to envMap_SH.cpp

// Calculates RGB Gradient (for color interpolation) and Luminance Hessian (for geometric error)
void calculateGradRGBAndHessianLumSHCoeffLM(
    const float3& x,
    const std::vector<ProbeSampleData>& samplingData,
    const std::vector<ProbeDirSample>& samplingDir,
    const int& basisIdx,
    GradSHCoeff& outGrad,
    float3x3& outHessian // Single 3x3 matrix for Luminance
)
{
    outGrad = { float3(0), float3(0), float3(0) };
    outHessian = float3x3::zeros();

    // Rec. 709 Luminance weights
    const float3 kLuma = float3(0.2126f, 0.7152f, 0.0722f);

    int samplingSize = samplingData.size();
    float Omega_i = (4.0f * M_PI) / (float)samplingSize;

    for (int sampleIdx = 0; sampleIdx < samplingSize; ++sampleIdx)
    {
        if (samplingData[sampleIdx].hitT < 0.0f) continue;

        float3 s = float3(samplingData[sampleIdx].s.x, samplingData[sampleIdx].s.y, samplingData[sampleIdx].s.z);
        float3 n = float3(samplingData[sampleIdx].n.x, samplingData[sampleIdx].n.y, samplingData[sampleIdx].n.z);

        // Input Radiance (RGB)
        float3 L = float3(samplingData[sampleIdx].Li.x, samplingData[sampleIdx].Li.y, samplingData[sampleIdx].Li.z);

        // Luminance Scalar for Hessian
        float L_lum = dot(L, kLuma);

        // --- Gradient (RGB) ---
        // We keep RGB gradients so runtime interpolation looks correct per-channel
        float3 gradOmega = gradientOmega(s, x, n, samplingSize);
        int numBasis = 9;
        float Ylm = SHBasisTable[numBasis * sampleIdx + basisIdx];
        float3 gradYlm = SHGradientTable[numBasis * sampleIdx + basisIdx];

        float3 q = s - x;
        float rInv = 1.0f / length(q);

        float3 contribGrad = gradOmega * Ylm - Omega_i * rInv * gradYlm;

        outGrad.r += (L.x * contribGrad);
        outGrad.g += (L.y * contribGrad);
        outGrad.b += (L.z * contribGrad);

        // --- Hessian (Luminance Only) ---
        // We use L_lum here. This saves 3x memory/compute vs storing RGB Hessians.
        float3x3 H_Omega = hessianOmega(s, x, n, samplingSize);
        float3x3 hessYlm = SHHessianTable[numBasis * sampleIdx + basisIdx];
        float3 w = samplingDir[sampleIdx].dir;
        float rInvSq = rInv * rInv;

        for (int j = 0; j < 3; ++j)
        {
            for (int k = 0; k < 3; ++k)
            {
                // Term 1: (∂²Ω / ∂xj∂xk) * Y
                float term1 = H_Omega[j][k] * Ylm;
                // Term 2: -1/r * [ (∂Ω/∂xj)(∂Y/∂wk) + (∂Ω/∂xk)(∂Y/∂wj) ]
                // Note: We access float3 components using array indexing [j] and [k]
                float term2 = -rInv * (gradOmega[j] * gradYlm[k] + gradOmega[k] * gradYlm[j]);
                // Term 3: Ω/r² * [ (ω)_j * (∂Y/∂wk) + (∂²Y / ∂wj∂wk) ]
                //float term3 = Omega_i * rInvSq * (w[j] * gradYlm[k] + hessYlm[j][k]);
                float term3 = Omega_i * rInvSq * hessYlm[j][k];

                // Accumulate Luminance Weighted Hessian
                outHessian[j][k] += L_lum * (term1 + term2 + term3);
            }
        }
    }
}

// Wrapper to loop over all 9 basis functions
void calculateSHCoeffsGradientsRGBAndHessiansLum(
    std::vector<GradSHCoeff>& gradOut,
    std::vector<float3x3>& hessianOut,
    const float3& gridPos,
    const std::vector<ProbeSampleData>& probeSamplingResults,
    const std::vector<ProbeDirSample>& samplingDir
)
{
    gradOut.resize(9);
    hessianOut.resize(9);

    for (int basisIdx = 0; basisIdx < 9; ++basisIdx)
    {
        calculateGradRGBAndHessianLumSHCoeffLM(
            gridPos,
            probeSamplingResults,
            samplingDir,
            basisIdx,
            gradOut[basisIdx],
            hessianOut[basisIdx]
        );
    }
}

//

// NEW: Compute SH coefficients for Radial Distance (r) and Squared Distance (r^2)
void calculateSHCoeffsRadialMoments(
    std::vector<float>& outMean,      // Output: 9 coeffs for E[r]
    std::vector<float>& outMeanSq,    // Output: 9 coeffs for E[r^2]
    const std::vector<ProbeSampleData>& probeSamplingResults,
    int numSamplePerProbe
)
{
    if (shOrder == -1 || SHBasisTable == nullptr) {
        // logError("Init SH table first");
        return;
    }

    int num_basis = 9;
    outMean.assign(num_basis, 0.0f);
    outMeanSq.assign(num_basis, 0.0f);

    float weight = 4.0f * float(M_PI) / float(numSamplePerProbe);

    for (int basisIdx = 0; basisIdx < num_basis; ++basisIdx)
    {
        double sumR = 0.0;
        double sumR2 = 0.0;

        for (int sampleIdx = 0; sampleIdx < numSamplePerProbe; ++sampleIdx)
        {
            float dist = probeSamplingResults[sampleIdx].hitT;

            // Handle Sky/Miss: Treat as "Far Away" (e.g., 100m)
            if (dist < 0.0f) dist = 10000.0f;

            float shBasis = SHBasisTable[num_basis * sampleIdx + basisIdx];

            sumR += (double)dist * shBasis * weight;
            sumR2 += (double)(dist * dist) * shBasis * weight;
        }

        outMean[basisIdx] = (float)sumR;
        outMeanSq[basisIdx] = (float)sumR2;
    }
}
