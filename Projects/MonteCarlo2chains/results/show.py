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

fileName = "results/Configurations/0confR.dat";
if(len(sys.argv)<2):
    confNum = 0;
else:
    confNum = int(sys.argv[1]);
#print fp.read(1);

#--------read and count chains and number of sites-------------
linesInBlock = 0;
prevLineIsEmpty = False;
emptyLine = "\n";

blockNum = 0;
lineNum=0;
N = [];
chain = [];
chain.append(Chain());

with open(fileName) as fp:
    for line in fp:
	if (line==emptyLine):
    	    if(prevLineIsEmpty==False):
        	N.append(lineNum - 1);
        	blockNum = blockNum + 1;
        	chain.append(Chain());
        	lineNum = 0;
    	    prevLineIsEmpty = True;
	else:
	    x, y, z = line.split();
	    chain[blockNum].pushCoordinates(float(x),float(y),float(z));
    	    lineNum = lineNum + 1;
    	    prevLineIsEmpty = False;
    
blockNum = blockNum-1;
#---------------------------end--------------------------------
print "Number of chains: ";
print blockNum;

fig = plt.figure()
ax = fig.gca(projection='3d');
ax.set_aspect('equal');

# Draw chain
ax.plot(chain[confNum].x, chain[confNum].y, chain[confNum].z, c='#006699');
ax.scatter(chain[confNum].x, chain[confNum].y, chain[confNum].z, c='#e60000', alpha=1);

# Make axes to be equal
X=np.asarray(chain[confNum].x);
Y=np.asarray(chain[confNum].y);
Z=np.asarray(chain[confNum].z);
max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max() / 2.0;
mid_x = (X.max()+X.min()) * 0.5;
mid_y = (Y.max()+Y.min()) * 0.5;
mid_z = (Z.max()+Z.min()) * 0.5;
ax.set_xlim(mid_x - max_range, mid_x + max_range);
ax.set_ylim(mid_y - max_range, mid_y + max_range);
ax.set_zlim(mid_z - max_range, mid_z + max_range);

plt.show();
