from ovito import dataset
import numpy as np
from ovito.vis import *
from ovito.io import import_file
from ovito.modifiers import *
import pylab as pl
from ovito.data import NearestNeighborFinder
from ovito.data import CutoffNeighborFinder
	




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



filexyz = 'bluePlots/current_excitations_T0.55_N10002_NVT_step_0.1LJ_startFrame500_tah200.xyz'
node=import_file(filexyz, multiple_frames = True,columns =[ "Particle Type", "Position.X", "Position.Y", "Position.Z"])
node.modifiers.append(ClusterAnalysisModifier(cutoff = 1.3, sort_by_size = True))
cutoff = 2
data = node.compute(500)
partType = 3
node.modifiers.append(SelectTypeModifier(types={partType}))
node.modifiers.append(InvertSelectionModifier())
node.modifiers.append(DeleteSelectedModifier())

node.modifiers.append(set_cell)
numberP = []
counterA = 0
counterB = 0
NNCount = []
#45.75777573740388 0.9769310226098933 for T= 0.6, strange, try other frames later
# 45.41906202723147 1.1800302571860817 for T= 0.5
# 46.02609188882587 0.7016449234259784 for T = 0.55
#29.99200799200799 0.7792207792207793
for frame in range(node.source.num_frames):
	#node.modifiers.append(randomParticleType)
	data = node.compute(frame)
	
	if frame%10 ==0:
		print('frame: ',frame)
	if data.particles.count>1:
		counterA_temp = countNeigh(data)
		NNCount.append(counterA_temp/data.particles.count)
		#print(frame,counterA_temp,data.particles.count)
		numberP.append(data.particles.count)
	elif data.particles.count==1:
		NNCount.append(0)
		numberP.append(1)
	

pl.hist(NNCount,normed = True,label = 'NN')
NNCount=np.array(NNCount)
numberP = np.array(numberP)
print(sum(NNCount==0),sum(NNCount==1),sum(NNCount>1),max(NNCount))
print(sum(numberP==0),sum(numberP==1),sum(numberP>1),max(numberP))
pl.hist(numberP,alpha=0.5,normed=True,label='numPart')
pl.legend(frameon=False)
np.save('NN_T055_exCount.npy',[NNCount,numberP])
pl.show()
