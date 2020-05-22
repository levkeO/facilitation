from ovito import dataset
import numpy as np
from ovito.vis import *
from ovito.io import import_file
from ovito.modifiers import *
import pylab as pl
from ovito.data import NearestNeighborFinder
from ovito.data import CutoffNeighborFinder
import pandas as pd
import sys
import singPartDist as sp





def set_cell(frame, data):
    with data.cell_:
        L = 19.25985167448
        data.cell_[:,0] = [L, 0., 0.]
        data.cell_[:,1] = [0., L, 0.]
        data.cell_[:,2] = [0., 0., L]
        #cell origin
        data.cell_[:,3] = [0,  0  ,  0]
        #set periodic boundary conditions
        data.cell_.pbc = (True, True, True)


# Prefetch the property array containing the particle type information:
def countNeigh(data):
	finder = CutoffNeighborFinder(cutoff, data)
	counter=0
	for index in range(data.particles.count):
		partCount =0
		for neigh in finder.find(index):
			counter+=1
	return counter


randSel = 2
def randomParticleType(frame,data):
        """
        Randomly reassignes particle types in the same ratio as before (for chosen particle type)
        Ovito modifier
        """
        randType = np.random.random_sample(size=(data.particles.count))
        randType= (randType<(randSel/data.particles.count)).astype(int)
        data.particles_.create_property('randProp', data=randType)



def pairdist(coords,p1,p2,frame1,frame2,L):
	"""
	Calculates the distance between two particles in a given frame range
	"""
	dist = coords[frame1:frame2,p2,:]-coords[frame1:frame2,p1,:]
	for frame in range(dist.shape[0]):
		dist[frame,:] = sp.periodic_boundary(dist[frame,:],L)
	return pl.sqrt(dist[:,0]**2+dist[:,1]**2+dist[:,2]**2)




excitations = pd.read_csv('excitation_results_T'+sys.argv[2]+'_tLJ01.csv')
for keys in excitations:
	if not keys == 'Unnamed: 0':
		node = import_file(sys.argv[1]+'T'+sys.argv[2]+'_N10002_NVT_step_0.1LJ_startFrame'+keys+'.xyz', multiple_frames=True,columns =["Particle Type", "Position.X", "Position.Y", "Position.Z"])
		partId = np.array(excitations[keys][0][1:-1].split(',')).astype(int)
		deltat = np.array(excitations[keys][1][1:-1].split(',')).astype(float)
		t0 = np.array(excitations[keys][2][1:-1].split(',')).astype(float)
		startFrame = (t0-deltat/2).astype(int)
		endFrame = (t0+deltat/2).astype(int)
		partId = partId[np.argsort(startFrame)]
		endFrame =  endFrame[np.argsort(startFrame)]

starts = []
snakes = {}
taken = []
rho=1.4
N = 10002
L  = (N/rho)**(1./3.)
numFrames = 1000
print(len(partId))
keys ='500'
coords=sp.readCoords(sys.argv[1]+'T'+sys.argv[2]+'_N10002_NVT_step_0.1LJ_startFrame'+keys+'.xyz',numFrames,N)
for index,particle in enumerate(partId):
	if not particle in taken:
		starts.append(particle)
		taken.append(particle)
		snakes[particle] = []
		followers = list(np.where(startFrame<endFrame[index])[0])
		for Id in followers:
			if not partId[Id] in taken:
				dist = pairdist(coords,particle,partId[Id],startFrame[Id],endFrame[index],L)
				if min(dist)<1.5:
					snakes[particle].append(partId[Id])
					starts.append(partId[Id])
					taken.append(partId[Id])

for item in snakes:
	if len(snakes[item])>0:
		print(item,snakes[item])
