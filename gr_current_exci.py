from ovito.io import *
from ovito.data import *
from ovito.modifiers import *
import numpy as np
import pylab as pl
import sys

node = import_file(sys.argv[1]+sys.argv[2],multiple_frames=True,columns =["Particle Type", "Position.X", "Position.Y", "Position.Z"])

def gr_partType(node,cut,bins,partType = 1,startFrame=0,endFrame=node.source.num_frames,step=1,randSel = None):
	"""
	Calculates the radial distribution function g(r) of one particle type averaged over all frames
	Only works with Ovito
	Args: 

		node (ovito data pipeline): all the data from file
		partType (int): index of the desired particle type
		cut (float): cutoff for g(r)
		bins(int): number of bins for g(r) histogram
	Returns:
		r (numpy array): r values of g(r)
		gr: radial distribution function g(r)for r values
		
	"""
	if randSel:
		node.modifiers.append(randomParticleType)
		node.modifiers.append(ExpressionSelectionModifier(expression ="randProp == 1"))
	
	else:
		node.modifiers.append(SelectTypeModifier(types={partType}))
	node.modifiers.append(InvertSelectionModifier())
	node.modifiers.append(DeleteSelectedModifier())
	node.modifiers.append(set_cell)
	modifier = CoordinationAnalysisModifier(cutoff = cut, number_of_bins = bins)
	node.modifiers.append(modifier)
	rdf = np.zeros((bins,2), float)
	counted = 0
	for frame in range(startFrame,endFrame,step):
		if frame%10 == 0:
			print (frame)
		# # Compute normalized bond vectors
		data = node.compute(frame)
		if not sum(data.particles['Particle Type'].array==3)>1:
			continue
		rdf+=data.tables['coordination-rdf'].xy()
		counted +=1
	rdf/=counted
	r = rdf[:,0]
	gr = rdf[:,1]
	print(counted)
	return r,gr
#node.modifiers.append(ExpressionSelectionModifier(expression ="ParticleType==1"))
#data = node.compute(1)
#numPart=data.particles.count
#print('number of particles: ',numPart)
#Assign random extra Type
#numType = np.count_nonzero(data.particles.selection) #why is this number not constant???
#print('number of excited particles',numType)


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



randSel = 850
def randomParticleType(frame,data):
	"""
	Randomly reassignes particle types in the same ratio as before (for chosen particle type)
	Ovito modifier
	"""
	randType = np.random.random_sample(size=(data.particles.count))
	randType= (randType<(randSel/data.particles.count)).astype(int)
	data.particles_.create_property('randProp', data=randType)
	

partType = 3

#r,gr = gr_partType(node,10,250,startFrame = 0,endFrame=500,randSel = randSel)
#pl.plot(r,gr, label ='random')
pl.rcParams.update({'font.size': 16})
node2 = import_file(sys.argv[1]+sys.argv[2],multiple_frames=True,columns =["Particle Type", "Position.X", "Position.Y", "Position.Z"])
data = node.compute(520)
r,gr = gr_partType(node2,10,250,startFrame=0,endFrame=1000,partType = partType)
pl.plot(r,gr,label = 'excited')
pl.legend(frameon=False)
#np.save(sys.argv[1]+"partType"+str(partType)+"_fast.npy", [r,gr])
pl.xlabel(r'$r$')
pl.ylabel(r'$g(r)$')

#node.modifiers.append(ExpressionSelectionModifier(expression ="ParticleType==1"))
#data = node.compute(1)
#numType = np.count_nonzero(data.particles.selection) #why is this number not constant???
#print('number of excited particles random',numType)
#pl.ylim([0,2])
pl.tight_layout()
#pl.savefig('gr_T06.pdf')
pl.show()
np.save('gr_current_'+sys.argv[2][:-4]+'.npy',[r,gr])
