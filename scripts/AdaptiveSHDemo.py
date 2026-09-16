from pathlib import WindowsPath, PosixPath
from falcor import *

def render_graph_AdaptiveSHDemo():
    g = RenderGraph('AdaptiveSHDemo')
    g.create_pass('AdaptiveSHDemo', 'AdaptiveSHDemo', {'runtimeGlass': True, 'glassRefraction': 1.0, 'liquidOpticalDepth': 0.5})
    g.create_pass('ToneMapper', 'ToneMapper', {'outputSize': 'Default', 'useSceneMetadata': True, 'exposureCompensation': 0.0, 'autoExposure': False, 'filmSpeed': 100.0, 'whiteBalance': False, 'whitePoint': 6500.0, 'operator': 'Reinhard', 'clamp': True, 'whiteMaxLuminance': 1.0, 'whiteScale': 11.199999809265137, 'fNumber': 1.0, 'shutter': 1.0, 'exposureMode': 'AperturePriority'})
    g.add_edge('AdaptiveSHDemo.output', 'ToneMapper.src')
    g.mark_output('AdaptiveSHDemo.output')
    g.mark_output('ToneMapper.dst')
    return g

AdaptiveSHDemo = render_graph_AdaptiveSHDemo()
try: m.addGraph(AdaptiveSHDemo)
except NameError: None
