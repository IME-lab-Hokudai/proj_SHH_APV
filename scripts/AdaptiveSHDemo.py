from pathlib import WindowsPath, PosixPath
from falcor import *

def render_graph_AdaptiveSHDemo():
    g = RenderGraph('AdaptiveSHDemo')
    g.create_pass('AdaptiveSHDemo', 'AdaptiveSHDemo', {})
    g.mark_output('AdaptiveSHDemo.output')
    return g

AdaptiveSHDemo = render_graph_AdaptiveSHDemo()
try: m.addGraph(AdaptiveSHDemo)
except NameError: None
