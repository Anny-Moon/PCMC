from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import sys

class Chain():
    def __init__(self):
	self.x =  [];
	self.y =  [];
	self.z =  [];
    def pushCoordinates(self, x, y, z):
	self.x.append(x);
	self.y.append(y);
	self.z.append(z);
    def getX(self):
	return self.x;
    def getY(self):
	return self.y;
    def getZ(self):
	return self.z;
class ChainFile():
    def __init__(self, fileName):
	#--------read and count chains and number of sites-------------
	linesInBlock = 0;
	prevLineIsEmpty = False;
	emptyLine = "\n";

	blockNum = 0;
	lineNum=0;
	self.N = [];
	self.chain = [];
	self.chain.append(Chain());

	with open(fileName) as fp:
	    for line in fp:
		if (line==emptyLine):
    		    if(prevLineIsEmpty==False):
        		self.N.append(lineNum - 1);
        		blockNum = blockNum + 1;
        		self.chain.append(Chain());
        		lineNum = 0;
    		    prevLineIsEmpty = True;
		else:
		    x, y, z = line.split();
		    self.chain[blockNum].pushCoordinates(float(x),float(y),float(z));
    		    lineNum = lineNum + 1;
    		    prevLineIsEmpty = False;
    
	self.numBlocks = blockNum-1;
	#---------------------------end--------------------------------
	print "Number of chains: ";
	print self.numBlocks;

    def plot(self, ax, confNum):
	ax.plot(self.chain[confNum].x, self.chain[confNum].y, self.chain[confNum].z, c='#006699');
	ax.scatter(self.chain[confNum].x, self.chain[confNum].y, self.chain[confNum].z, c='#e60000', alpha=1);
    
    def getChain(self, chainNum):
	return self.chain[chainNum];
	
class EqualAxes():
    def __init__(self):
	self.X = [];
	self.Y = [];
	self.Z = [];
    def push(self,X,Y,Z):
	self.X.append(X);
	self.Y.append(Y);
	self.Z.append(Z);
    def set(self,ax):
	X=np.asarray(self.X);
	Y=np.asarray(self.Y);
	Z=np.asarray(self.Z);
	max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max() / 2.0;
	mid_x = (X.max()+X.min()) * 0.5;
	mid_y = (Y.max()+Y.min()) * 0.5;
	mid_z = (Z.max()+Z.min()) * 0.5;
	ax.set_xlim(mid_x - max_range, mid_x + max_range);
	ax.set_ylim(mid_y - max_range, mid_y + max_range);
	ax.set_zlim(mid_z - max_range, mid_z + max_range);


fileName = "Configurations/0confR.dat";
chain = ChainFile(fileName);
	
if(len(sys.argv)<2):
    confNum = 0;
else:
    confNum = int(sys.argv[1]);

fig = plt.figure()
ax = fig.gca(projection='3d');
ax.set_aspect('equal');
chain.plot(ax,confNum);
eqAx = EqualAxes();
eqAx.push(chain.getChain(confNum).getX(),chain.getChain(confNum).getY(),chain.getChain(confNum).getZ());
eqAx.set(ax);

plt.show();
