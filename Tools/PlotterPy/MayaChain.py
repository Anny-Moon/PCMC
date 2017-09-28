from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import random
import Color
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from mayavi import mlab

class FixZorderCollection(Line3DCollection):
    _zorder = 1000

    @property
    def zorder(self):
    	return self._zorder

    @zorder.setter
    def zorder(self, value):
    	pass

class Chain():
    def __init__(self):
	self.x =  [];
	self.y =  [];
	self.z =  [];

    def pushCoordinates(self, x, y, z):
	self.x.append(x);
	self.y.append(y);
	self.z.append(z);
    
    def plot(self, fig, sizeDot=None, sizeLine=None, colorDot=None, colorLine=None):
	
	mlab.points3d(self.x,self.y,self.z, color = colorDot, scale_factor=sizeDot, figure = fig);
	mlab.plot3d(self.x,self.y,self.z, tube_radius=sizeLine, color = colorLine, figure = fig);
	
	
    def happyPlot(self, ax, axMaxRange=None):
	if(axMaxRange == None):
	    axMaxRange = 20;
	sizeLine = 30/axMaxRange;
	if(sizeLine<1):
	    sizeLine=1;
	sizeDot = 800/axMaxRange;
	#sizeDot = 1000/axMaxRange;
	
	colorLine = "#000000";

	N=len(self.x);
	for i in range (0, N-1):
	    r = lambda: random.randint(0,255);
	    colorDot = '#%02X%02X%02X' % (r(),r(),r());
	    
	    line = [(self.x[i],self.y[i],self.z[i]),
		(self.x[i+1],self.y[i+1],self.z[i+1])];
	    (x, y, z) = zip(*line);
	    ax.scatter(x,y,z, c = colorDot, alpha = 1, s=sizeDot);
	    ax.plot(x,y,z, c = colorLine, alpha = 1, lw=sizeLine);
	
    def plotOld(self, ax, colorDot = None, colorLine = None):
	if(colorDot == None):
	    colorDot = '#e60000';
	if(colorLine == None):
	    colorLine = '#006699';
	ax.plot(self.x, self.y, self.z, c=colorLine, zorder=1);
	ax.scatter(self.x, self.y, self.z, c=colorDot, alpha=1, s = 50, zorder=2);

    def getX(self):
	return self.x;
    def getY(self):
	return self.y;
    def getZ(self):
	return self.z;
    def getNumSites(self):
	return len(self.x);