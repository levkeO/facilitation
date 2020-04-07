"""
function to detect excitations and caclulate the distribution of Delta t (and later a)
"""

import sys
import glob
path='/mnt/e/Simulations/Facilitation/facilitation/'
sys.path.append(path)	
import singPartDist as sp
from numba import njit, config, __version__
from numba.extending import overload
import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit
import pylab as pl

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

def detect(fileName, tah, a, N,numFrames, L):
	"""
	Detects the number of excitations in the given file
	Args:
		fileName (string): name of xyz file
		tah (int): timescale over which to search for excitations
		a (float): lenght scale as threshold for excitations
		N (int): number of particles
		numFrames (int): number of frames
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
	block = sp.readCoords(fileName,numFrames,N)
	for particle in range(N):
		p_coord = block[:,particle,:]
		p_deltat = []
		t0_temp = []
		cond1 =sp.averageDistPos(p_coord,0,tah,-1-tah,-1,0,L) 
		if cond1 > (2*a)**2:
			for frame in range(tah,numFrames-tah):
				cond2 = sp.averageDistPos(p_coord,frame-tah,frame,frame,frame+tah,frame,L)
				if cond2 > (2*a)**2:
					p = sp.singlepath(p_coord,frame,L)
					t = np.arange(p.diffs.shape[0])
					w = np.append(0,np.where(p.diffs[frame-tah:frame+tah]-a<0)[0])
					last = max(w)+1 # last time the distance drops below a/2	
					first = min(w[w>0]) # last time in other direction the  distance drops below a/2
					p_deltat.append(last-first) # append deltat to particle array
					t0_temp.append(frame)
				if len(p_deltat)>0:
					p_t0  = t0_temp[np.argmin(p_deltat)]
					exPart.append(particle)
					t_fit = t[p_t0-tah:p_t0+tah]
					p_fit = p.diffs[p_t0-tah:p_t0+tah]
					spl =UnivariateSpline(t_fit,p_fit,s=0.8)
					pl.figure()
					pl.plot(t_fit,p_fit)
					pl.plot(t_fit,spl(t_fit))
					pl.title(particle)
					try:
						fit_param, covar = curve_fit(mytanh,t_fit,spl(t_fit),[tah/20,0.3,p_t0,0.5], bounds = ([0,-np.inf,0,-np.inf], [np.inf,np.inf,np.inf,np.inf]))
					except:
						print('Particle ',particle,'failed!')
						break
					deltat.append(fit_param[0]*2)
					start.append(fit_max[0][0])		
					t0.append(fit_param[2])
	pl.show()
	
	return exPart, deltat, start, t0
					




fileName = glob.glob('../Data/T0.5/*')
print(fileName[0])
rho=1.4
N = 10002
L  = (N/rho)**(1./3.)
a,b,c,d = detect(fileName[0],200,0.15,N,1000,L)
#print(min(b),len(b))
pl.show()


