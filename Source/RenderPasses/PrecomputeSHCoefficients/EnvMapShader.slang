#include "Utils/Math/MathConstants.slangh"
import Scene.Shading;

struct VsOut
{
    float2 texC : TEXCOORD;
#ifndef _VIEWPORT_MASK
    float4 posH : SV_POSITION;
#else
    float4 posH : POSITION;
#endif
};

cbuffer PerFrameCB
{
    float4 shCoeffs[9];
    bool showReconstructedEnvMap;
};

// Evaluate 9 SH basis functions at direction n
float4 evaluateIrradiance(float3 n)
{
    float x = n.x;
    float y = n.y;
    float z = n.z;

    // SH basis function indexing and formulas:
    // Index | l (band) | m (order) | SH function          | Formula (Cartesian coords)
    // -------------------------------------------------------------------------------------
    // 0     | 0        |  0        | Y_0^0                | 0.282095
    // 1     | 1        | -1        | Y_1^-1               | 0.488603 * y       // ~ sqrt(3/(4π)) * y/r * sin(φ)
    // 2     | 1        |  0        | Y_1^0                | 0.488603 * z       // ~ sqrt(3/(4π)) * z/r
    // 3     | 1        |  1        | Y_1^1                | 0.488603 * x       // ~ sqrt(3/(4π)) * x/r * cos(φ)
    // 4     | 2        | -2        | Y_2^-2               | 1.092548 * x * y   // ~ sqrt(15/(4π)) * xy/r² * sin(2φ)
    // 5     | 2        | -1        | Y_2^-1               | 1.092548 * y * z   // ~ sqrt(15/(4π)) * yz/r² * sin(φ)
    // 6     | 2        |  0        | Y_2^0                | 0.315392 * (-x*x - y*y + 2*z*z) // ~ 0.315392 * (3z² - 1)
    // 7     | 2        |  1        | Y_2^1                | 1.092548 * x * z   // ~ sqrt(15/(4π)) * xz/r² * cos(φ)
    // 8     | 2        |  2        | Y_2^2                | 0.546274 * (x*x - y*y)  // ~ 0.546274 * (x² - y²)/r² * cos(2φ)

    // Band 0
    float Y0 = 0.282095;

    // Band 1
    float Y1 = 0.488603 * y;
    float Y2 = 0.488603 * z;
    float Y3 = 0.488603 * x;

    // Band 2
    float Y4 = 1.092548 * x * y;
    float Y5 = 1.092548 * y * z;
    float Y6 = 0.315392 * (3.0 * z * z - 1.0);
    float Y7 = 1.092548 * x * z;
    float Y8 = 0.546274 * (x * x - y * y);

    // Dot product: coeffs · basis
    float4 result = 
        shCoeffs[0] * Y0 +
        shCoeffs[1] * Y1 +
        shCoeffs[2] * Y2 +
        shCoeffs[3] * Y3 +
        shCoeffs[4] * Y4 +
        shCoeffs[5] * Y5 +
        shCoeffs[6] * Y6 +
        shCoeffs[7] * Y7 +
        shCoeffs[8] * Y8;

    return max(result, float4(0.0)); // Clamp to avoid negative lighting
}

VsOut vsMain(float4 posS: POSITION, float2 texC: TEXCOORD)
{
    VsOut vOut;
    vOut.texC = texC;
    vOut.posH = posS;

    return vOut;
}

float4 main(VsOut vsOut) : SV_TARGET
{

    // Convert screen-space UV to NDC
    float2 ndc = vsOut.texC * 2.0 - 1.0;
    float4 clipPos = float4(ndc.x, -ndc.y, 1.0, 1.0); // z = 1 for far plane 

    // Get camera parameters
    float4x4 gInvViewProj = gScene.camera.data.invViewProj;
    float4x4 viewMat = gScene.camera.data.viewMat;

    // Reconstruct world-space position on far plane
    float4 worldPosH = mul(gInvViewProj, clipPos);
    float3 worldPos = worldPosH.xyz / worldPosH.w;

    // Step 3: Reconstruct world-space ray direction
    float3 dir = normalize(worldPos - gScene.camera.getPosition());

    if (showReconstructedEnvMap) {
        return float4(evaluateIrradiance(dir).rgb, 1.0f);
    }

    float2 uv = gScene.envMap.worldToUv(dir);
    float3 color = gScene.envMap.eval(uv);
    
    return float4(color, 1.0f);
}
