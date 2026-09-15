add_plugin(BakeLightMapXAtlas)

# The baker imports the original PathTracer shader modules at runtime.
add_dependencies(BakeLightMapXAtlas PathTracer)

target_sources(BakeLightMapXAtlas PRIVATE
    BakeLightMapXAtlas.cpp
    AtlasSeams.cpp
    BakeLightMapXAtlas.h
    xatlas.h
    xatlas.cpp
    XAtlasUVRaster.slang
	BakeDataStructures.slang
	ExtractTexelsSingle.cs.slang
	LightmapBakingSingle.rt.slang
	NormalizeLightmapSingle.cs.slang
    BlurLightmap.cs.slang
    DilateLightmap.cs.slang
    AtlasSeams.cs.slang
)
target_copy_shaders(BakeLightMapXAtlas RenderPasses/BakeLightMapXAtlas)
target_source_group(BakeLightMapXAtlas "RenderPasses")
