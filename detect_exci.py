"""
function to detect excitations and caclulate the distribution of Delta t (and later a)
"""

import sys
import glob
path='/mnt/e/Simulations/Facilitation/facilitation/'
sys.path.append(path)	
import singPartDist as sp

def tanh(x,c1,t0,a,c2):
	"""
	hyberbolic tangens with for parameters
	
	Args:
		x (list or numpy array or float): input parameter for tanh
		c1 (float): multiplicative parameter  
		t0 (float): parameter to shift the centre
		a (float):  stretche or squeeze tanh in x-direction
		c2 (float): additive parameter, shift tanh up and down	
	Return:
		value of c1* tanh((x-t0)/a)+c2
	"""
	return c1*np.tanh((x-t0)/a)+c2

def detect(fileName, tah, a, N, L)
	"""
	Detects the number of excitations in the given file
	Args:
		fileName (string): name of xyz file
		tah (int): timescale over which to search for excitations
		a (float): lenght scale as threshold for excitations
		N (int): number of particles
		L (float):  length of box
	Return:
		array of all Delta t
		startT
		expart ID
		t0s
	"""
	start = []
	deltat = []
	exPart = []
	t0 = []
	block = sp.readCoords(fileName,numFrames,numPart)
	for particle in range(N):
		p_coord = block[:,particle,:]
		p_start = []
		p_deltat = []
		p_t0 = []
		cond1 =sp.averageDistPos(p_coord,0,tah,-1-tah,-1,0,L) 
		if cond > (2*a)**2:
			for frame in range(tah,numFrames-tah):
				cond2 = sp.averageDistPos(p_coord,frame-tah,frame,frame,frame+tah,frame,L)















