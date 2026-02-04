from pathlib import WindowsPath, PosixPath
from falcor import *

def render_graph_BakeLightMap():
    g = RenderGraph('BakeLightMap')
    g.create_pass('BakeLightMap', 'BakeLightMap', {})
    g.mark_output('BakeLightMap.output')
    return g

BakeLightMap = render_graph_BakeLightMap()
try: m.addGraph(BakeLightMap)
except NameError: None
