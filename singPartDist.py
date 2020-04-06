"""Module to calculate distance measures of single particle trajectories

"""

import pylab as pl
import numpy as np
import sys
from numba import njit, config, __version__
from numba.extending import overload


def readCoords(filexyz, numFrames, numPart):
	"""
	Reads data from an xyz file
	Args:
		filexyz(string): name of the xyz file to read
		numFrames (int): number of frames in file
		numPart (int): number of particles
	Return:
		allCoords (list of list) for each frame a list of all 
		particles consisting of list of all three coordinates
		for each particle (x,y,z)
	
	"""

	frame = -1
	allCoords = np.zeros((numFrames,numPart,3))
	with open(filexyz, 'r') as readFile:
		for line in readFile:
			splitL = line.split()
			if len(splitL) ==1:
				frame +=1
				if frame == numFrames:
					break
				particleCounter = 0
			elif not splitL[0] == 'Atoms.':
				allCoords[frame][particleCounter,0] =splitL[1]
				allCoords[frame][particleCounter,1] =splitL[2]
				allCoords[frame][particleCounter,2] =splitL[3]
				particleCounter+=1
	return allCoords


@njit
def periodic_boundary(xyzArray,L):
        """     
        Makes sure that the given coordinates are inside the box of Length L (between -L/2 and L/2) 
        And applies periodic boundary conditions
	Origin is in the middle of the box, could be an option to put it at the left later
        
        Args:
        	xyzArray (array with 3 entries): array with one set of coordinates (3D)
        	L(int):, length of box, square box assumed

        Return:
        	an xyz-coordinate array inside the box
	Examples:
		>>> periodic_boundary([1,3,6],10)
		[1, 3, -4]
        """
        for dim in range(3):
                if xyzArray[dim]>L/2 : xyzArray[dim]-=L
                if xyzArray[dim]<-L/2 : xyzArray[dim]+=L
        return xyzArray


@njit
def squareDist(coords, frame1, frame2,L):
	"""
	Distance between two frames of one particle taking care of boundary 
	conditions
	Args:
		coords (list of list): coordinates of one particle for several frames in  3d
		frame1 (int): index of first frame to calculate distance 
		frame2 (int): index of second frame to calculate distance 
        	L(int):, length of box, square box assumed
	Return:
		squared distance (float)
	Example:
		dsjkldsyjdsyjkl
	"""
	dist = coords[frame2,:]-coords[frame1,:]
	dist = periodic_boundary(dist,L)
	return dist[0]**2 + dist[1]**2 + dist[2]**2


@njit
def averageDistPos(coords, start1,end1,start2,end2,reference,L):
	# not sure if this is a sensible setup, probably not
	dist1 = coords[start1:end1,:]- coords[reference,:]
	dist2 = coords[start2:end2,:]- coords[reference,:]
	for t in range(len(dist1)):
		dist1[t] = periodic_boundary(dist1[t],L)
		dist2[t] = periodic_boundary(dist2[t],L)
	average1 = dist1.mean(axis=0)
	average2 = dist2.mean(axis=0)
	distance = average2-average1
	distance = periodic_boundary(distance,L)
	return distance[0]**2 +distance[1]**2 + distance[2]**2



		