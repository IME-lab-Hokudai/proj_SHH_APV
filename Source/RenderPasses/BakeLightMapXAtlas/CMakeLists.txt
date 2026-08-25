add_plugin(BakeLightMapXAtlas)

target_sources(BakeLightMapXAtlas PRIVATE
    BakeLightMapXAtlas.cpp
    BakeLightMapXAtlas.h
    xatlas.h
    xatlas.cpp
    XAtlasUVRaster.slang
	BakeDataStructures.slang
	ExtractTexelsSingle.cs.slang
	LightmapBakingSingle.rt.slang
	NormalizeLightmapSingle.cs.slang
)
target_copy_shaders(BakeLightMapXAtlas RenderPasses/BakeLightMapXAtlas)
target_source_group(BakeLightMapXAtlas "RenderPasses")
