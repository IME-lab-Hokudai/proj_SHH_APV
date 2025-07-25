#include "Utils/Math/MathConstants.slangh"
import Scene.Shading;

float4 main(float2 texCoord: TEXCOORD) : SV_TARGET
{

    // Convert screen-space UV to NDC
    float2 ndc = texCoord * 2.0 - 1.0;
    float4 clipPos = float4(ndc.x, -ndc.y, 1.0, 1.0); // z = 1 for far plane

    // Get camera parameters
    float4x4 gInvViewProj = gScene.camera.data.invViewProj;
    float4x4 viewMat = gScene.camera.data.viewMat;

    // Reconstruct world-space position on far plane
    float4 worldPosH = mul(gInvViewProj, clipPos);
    float3 worldPos = worldPosH.xyz / worldPosH.w;

    // Step 3: Reconstruct world-space ray direction
    float3 dir = normalize(worldPos - gScene.camera.getPosition());

    float phi = atan2(dir.z, dir.x);
    float theta = acos(clamp(dir.y, -1.0, 1.0));
    phi = (phi < 0.0) ? phi + 2.0 * M_PI : phi;

    float2 uv = float2(phi / (2.0 * M_PI), theta / M_PI);

    float3 color = gScene.envMap.eval(uv);

    return float4(color, 1.0f);
}
