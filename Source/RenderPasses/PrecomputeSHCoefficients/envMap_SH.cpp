#include "sphericalHarmonics.h"
#include "envMap_SH.h"

float*  dOmega;
float*  SHBasisTable;
float3* SHGradientTable;
float3x3* SHHessianTable;
int shOrder = -1;

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

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

void calculateGradSHCoeffLM(
    const float3& x,
    const std::vector<ProbeSampleData>& samplingData,
    const std::vector<ProbeDirSample>& samplingDir,
    const int& basisIdx,
    GradSHCoeff& outGrad
)
{

    outGrad = { float3(0), float3(0), float3(0) };

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

        // q = s - x
        float3 q = s - x;

        float rInv = 1.0f / length(q);
        float cosXi = -(dot(n, q)) * rInv;

        float3 gradOmega = gradientOmega(q, n, rInv, cosXi, Omega_i);

        // SH basis and derivative in direction space
        int numBasis = 9;
        float Ylm = SHBasisTable[numBasis * sampleIdx + basisIdx];
        float3 gradYlm = SHGradientTable[numBasis * sampleIdx + basisIdx];

        // Contribution of this sample to ∇f_l^m
        // ∇_x Ω_i Y_l^m (ω_i )-Ω_i/r_i  ∇_(ω_i ) Y_l^m (ω_i)
        float3 contrib = gradOmega * Ylm - Omega_i * rInv * gradYlm;
        outGrad.r += (L.r * contrib);
        outGrad.g += (L.g * contrib);
        outGrad.b += (L.b * contrib);
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

float3 gradientOmega(const float3& q, const float3& n, float rInv, float cosXi, float factor)
{
    // Numerical guard for cosXi
    const float eps = 1e-6f;
    if (fabs(cosXi) < eps)
        cosXi = (cosXi < 0 ? -eps : eps);

    float invCosXi = 1.0f / cosXi;
    float rInv2 = rInv * rInv;

    // Analytically simplified from: factor * ((n * r + 3 * q * cosXi) / (r^2 * cosXi))
    // to: factor * (n / (r * cosXi) + 3 * q / r^2)
    return factor * (n * (rInv * invCosXi) + q * (3.0f * rInv2));
}

float3x3 hessianOmega(const float3& q, const float3& n, float rInv, float cosXi, float factor)
{
    float qx = q.x, qy = q.y, qz = q.z;

    // Numerical guard for cosXi
    const float eps = 1e-6f;
    if (fabs(cosXi) < eps)
        cosXi = (cosXi < 0 ? -eps : eps);

    // Precompute inverse powers (Multiplication is vastly faster than division)
    float invCosXi = 1.0f / cosXi;
    float rInv2 = rInv * rInv;
    float rInv3 = rInv2 * rInv;
    float rInv4 = rInv2 * rInv2;

    // Pure second derivatives
    float term_xx = 6.0f * n.x * qx * rInv3 * invCosXi - 3.0f * rInv2 + 15.0f * qx * qx * rInv4;
    float term_yy = 6.0f * n.y * qy * rInv3 * invCosXi - 3.0f * rInv2 + 15.0f * qy * qy * rInv4;
    float term_zz = 6.0f * n.z * qz * rInv3 * invCosXi - 3.0f * rInv2 + 15.0f * qz * qz * rInv4;

    // Mixed derivatives
    float term_xy = 3.0f * (n.x * qy + qx * n.y) * rInv3 * invCosXi + 15.0f * qx * qy * rInv4;
    float term_xz = 3.0f * (n.x * qz + qx * n.z) * rInv3 * invCosXi + 15.0f * qx * qz * rInv4;
    float term_yz = 3.0f * (n.y * qz + qy * n.z) * rInv3 * invCosXi + 15.0f * qy * qz * rInv4;

    // Assemble symmetric Hessian
    float3x3 H = float3x3::zeros();
    H[0][0] = factor * term_xx;
    H[0][1] = factor * term_xy;
    H[0][2] = factor * term_xz;

    H[1][0] = factor * term_xy;
    H[1][1] = factor * term_yy;
    H[1][2] = factor * term_yz;

    H[2][0] = factor * term_xz;
    H[2][1] = factor * term_yz;
    H[2][2] = factor * term_zz;

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

void calculateChannelRGradAndHessianSHCoeffLM(const float3& x, const std::vector<ProbeSampleData>& samplingData, const std::vector<ProbeDirSample>& samplingDir, const int& basisIdx, float3& outGrad, float3x3& outHessian)
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
        // q = s - x
        float3 q = s - x;

        // r = ||q||
        float rInv = 1.0f / length(q);
        float cosXi = -(dot(n, q)) * rInv;
        // Gradient ∂_x Ω_i
        float3 gradOmega = gradientOmega(q, n, rInv, cosXi, Omega_i);
        float3x3 H_Omega = hessianOmega(q, n, rInv, cosXi, Omega_i);
        // SH basis and derivative in direction space
        int numBasis = 9;
        float Ylm = SHBasisTable[numBasis * sampleIdx + basisIdx];
        float3 gradYlm = SHGradientTable[numBasis * sampleIdx + basisIdx];

        // Contribution of this sample to ∇f_l^m
        // ∇_x Ω_i Y_l^m(ω_i) - (Ω_i / r_i) ∇_{ω_i} Y_l^m(ω_i)
        float3 gradContrib = gradOmega * Ylm - Omega_i * rInv * gradYlm;
        outGrad += (L.r * gradContrib);

        float3x3 hessYlm = SHHessianTable[numBasis * sampleIdx + basisIdx];

        // Pre-calculate shared values for the inner loop
        float rInvSq = rInv * rInv;

        // Accumulate Hessian of f_l^m
        // Utilizing Schwarz's symmetry (H[j][k] == H[k][j]) to halve loop cycles
        for (int j = 0; j < 3; ++j)
        {
            for (int k = j; k < 3; ++k) // Notice k starts at j now
            {
                // Term 1: (∂²Ω / ∂x_j ∂x_k) * Y_l^m
                float term1 = H_Omega[j][k] * Ylm;

                // Term 2: -1/r * [ (∂Ω/∂x_j)(∂Y/∂ω_k) + (∂Ω/∂x_k)(∂Y/∂ω_j) ]
                float term2 = -rInv * (gradOmega[j] * gradYlm[k] + gradOmega[k] * gradYlm[j]);

                // Term 3: (Ω / r²) * (∂²Y / ∂ω_j ∂ω_k)
                float term3 = Omega_i * rInvSq * hessYlm[j][k];

                float hessContrib = term1 + term2 + term3;

                // Assign to the upper triangle
                outHessian[j][k] += L.r * hessContrib;

                // Mirror to the lower triangle (avoiding double-counting the diagonal)
                if (j != k)
                {
                    outHessian[k][j] += L.r * hessContrib;
                }
            }
        }
    }
}

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

        // q = s - x
        float3 q = s - x;

        float rInv = 1.0f / length(q);
        float cosXi = -(dot(n, q)) * rInv;

        float3 gradOmega = gradientOmega(q, n, rInv, cosXi, Omega_i);
        float3x3 H_Omega = hessianOmega(q, n, rInv, cosXi, Omega_i);

        int numBasis = 9;
        float Ylm = SHBasisTable[numBasis * sampleIdx + basisIdx];
        float3 gradYlm = SHGradientTable[numBasis * sampleIdx + basisIdx];

        float3 contribGrad = gradOmega * Ylm - Omega_i * rInv * gradYlm;

        outGrad.r += (L.x * contribGrad);
        outGrad.g += (L.y * contribGrad);
        outGrad.b += (L.z * contribGrad);

        // --- Hessian (Luminance Only) ---
        float3x3 hessYlm = SHHessianTable[numBasis * sampleIdx + basisIdx];
        float rInvSq = rInv * rInv;

        // Utilizing Schwarz's symmetry (H[j][k] == H[k][j]) to halve loop cycles
        for (int j = 0; j < 3; ++j)
        {
            for (int k = j; k < 3; ++k) // Notice k starts at j now
            {
                // Term 1: (∂²Ω / ∂x_j ∂x_k) * Y_l^m
                float term1 = H_Omega[j][k] * Ylm;

                // Term 2: -1/r * [ (∂Ω/∂x_j)(∂Y/∂ω_k) + (∂Ω/∂x_k)(∂Y/∂ω_j) ]
                float term2 = -rInv * (gradOmega[j] * gradYlm[k] + gradOmega[k] * gradYlm[j]);

                // Term 3: (Ω / r²) * (∂²Y / ∂ω_j ∂ω_k)
                float term3 = Omega_i * rInvSq * hessYlm[j][k];

                float hessContrib = term1 + term2 + term3;

                // Assign to the upper triangle
                outHessian[j][k] += L_lum * hessContrib;

                // Mirror to the lower triangle (avoiding double-counting the diagonal)
                if (j != k)
                {
                    outHessian[k][j] += L_lum * hessContrib;
                }
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
