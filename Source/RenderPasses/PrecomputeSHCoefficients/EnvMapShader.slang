#include "Utils/Math/MathConstants.slangh"
import Scene.Shading;

float4 main(float2 texCoord: TEXCOORD) : SV_TARGET
{
    // Convert 2D texCoord to direction for env map lookup
    float phi = texCoord.x * 2.0 * M_PI;
    float theta = texCoord.y * M_PI;

    float3 dir = float3(
        sin(theta) * cos(phi),
        cos(theta),
        sin(theta) * sin(phi)
    );

    float3 color = gScene.envMap.eval(dir);

    return float4(color, 1.0f);
}
