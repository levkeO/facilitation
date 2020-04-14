import sys
sys.path.append('/usr/local/lib/python3.5/dist-packages')
from ovito import dataset
import numpy as np
from ovito.vis import *
from ovito.io import import_file
from ovito.modifiers import *
import pylab as pl
from ovito.data import NearestNeighborFinder
from ovito.data import CutoffNeighborFinder
filexyz = '../../KA21/T0.5/excitations_T0.5_N10002_NVT_step_0.1LJ_startFrame450_tah200.xyz'
import_file(filexyz, multiple_frames = True,columns =[ "Particle Type", "Position.X", "Position.Y", "Position.Z"])
node.modifiers.append(ClusterAnalysisModifier(cutoff = 1.3, sort_by_size = True))

