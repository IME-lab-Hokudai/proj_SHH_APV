from pathlib import WindowsPath, PosixPath
from falcor import *

def render_graph_BakeLightMapXAtlas():
    g = RenderGraph('BakeLightMapXAtlas')
    g.create_pass('BakeLightMapXAtlas', 'BakeLightMapXAtlas', {})
    g.mark_output('BakeLightMapXAtlas.output')
    return g

BakeLightMapXAtlas = render_graph_BakeLightMapXAtlas()
try: m.addGraph(BakeLightMapXAtlas)
except NameError: None
