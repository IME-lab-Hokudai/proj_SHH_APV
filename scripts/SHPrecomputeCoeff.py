from pathlib import WindowsPath, PosixPath
from falcor import *

def render_graph_SHPrecompute():
    g = RenderGraph('SHPrecompute')
    g.create_pass('PrecomputeSHCoefficients', 'PrecomputeSHCoefficients', {})
    g.mark_output('PrecomputeSHCoefficients.output')
    return g

SHPrecompute = render_graph_SHPrecompute()
try: m.addGraph(SHPrecompute)
except NameError: None
