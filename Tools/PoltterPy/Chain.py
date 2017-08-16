from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import random

class Chain():
    def __init__(self):
	self.x =  [];
	self.y =  [];
	self.z =  [];

    def pushCoordinates(self, x, y, z):
	self.x.append(x);
	self.y.append(y);
	self.z.append(z);
    
    def plot(self, ax, axMaxRange=None, colorDot=None, colorLine=None):
	
	if(axMaxRange == None):
	    axMaxRange = 20;
	sizeLine = 30/axMaxRange;
	if(sizeLine<2):
	    sizeLine=2;
	#sizeDot = 800/axMaxRange;
	sizeDot = 1000/axMaxRange;
	
	if(colorDot == None):
	    r = lambda: random.randint(0,255);
	    r = lambda: random.randint(0,255);
	    colorDot = '#%02X%02X%02X' % (r(),r(),r());
	
	if(colorLine == None):
	    r = lambda: random.randint(0,255);
	    colorLine = '#%02X%02X%02X' % (r(),r(),r());

	N=len(self.x);
	for i in range (0, N-1):
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
