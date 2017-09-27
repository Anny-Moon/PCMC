import numpy as np

class EqualAxes():
    def __init__(self, ax):
	self.axes = ax;
	self.X = [];
	self.Y = [];
	self.Z = [];
	
    def push(self,X,Y,Z):
	self.X+=X;
	self.Y+=Y;
	self.Z+=Z;
	
    def findMaxRange(self):
	X=np.asarray(self.X);
	Y=np.asarray(self.Y);
	Z=np.asarray(self.Z);
	maxRange = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max();
	return maxRange;
	
    def set(self):
	X=np.asarray(self.X);
	Y=np.asarray(self.Y);
	Z=np.asarray(self.Z);
	max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max() / 2.0;
	mid_x = (X.max()+X.min()) * 0.5;
	mid_y = (Y.max()+Y.min()) * 0.5;
	mid_z = (Z.max()+Z.min()) * 0.5;
	self.axes.set_xlim(mid_x - max_range, mid_x + max_range);
	self.axes.set_ylim(mid_y - max_range, mid_y + max_range);
	self.axes.set_zlim(mid_z - max_range, mid_z + max_range);
    
    def getAxes(self):
	return self.axes;