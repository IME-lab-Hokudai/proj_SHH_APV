add_plugin(BakeLightMap)

target_sources(BakeLightMap PRIVATE
    BakeLightMap.cpp
    BakeLightMap.h
	LightmapBaking.rt.slang
	UVPass.slang
	BakeDataStructures.slang
	ExtractTexels.cs.slang
	NormalizeLightmap.cs.slang
	ApplyLightmap.slang
)

target_copy_shaders(BakeLightMap RenderPasses/BakeLightMap)
target_source_group(BakeLightMap "RenderPasses")
