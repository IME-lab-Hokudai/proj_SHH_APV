from pathlib import WindowsPath, PosixPath
from falcor import *

def render_graph_BakeLightMapSingle():
    g = RenderGraph('BakeLightMapSingle')
    g.create_pass('BakeLightMapSingle', 'BakeLightMapSingle', {})
    g.mark_output('BakeLightMapSingle.output')
    return g

BakeLightMapSingle = render_graph_BakeLightMapSingle()
try: m.addGraph(BakeLightMapSingle)
except NameError: None
