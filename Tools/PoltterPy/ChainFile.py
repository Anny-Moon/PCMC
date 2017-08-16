import Chain

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
	self.chain.append(Chain.Chain());

	with open(fileName) as fp:
	    for line in fp:
		if (line==emptyLine):
    		    if(prevLineIsEmpty==False):
        		self.N.append(lineNum);
        		blockNum = blockNum + 1;
        		self.chain.append(Chain.Chain());
        		lineNum = 0;
    		    prevLineIsEmpty = True;
		else:
		    x, y, z = line.split();
		    self.chain[blockNum].pushCoordinates(float(x),float(y),float(z));
    		    lineNum = lineNum + 1;
    		    prevLineIsEmpty = False;
    
	self.numBlocks = blockNum-1;
	#---------------------------end--------------------------------
	print('Number of chains: %i.' % self.numBlocks);
	
    def plot(self, confNum, ax, axMaxRange=None, colorDot=None, colorLine=None):
	self.chain[confNum].plot(ax, axMaxRange, colorDot, colorLine);
	
    def plotOld(self, confNum, ax, colorDot=None, colorLine=None):
	self.chain[confNum].plotOld(ax, colorDot, colorLine);
	
    def getChain(self, chainNum):
	return self.chain[chainNum];
    def getN(self, chainNum):
	return self.N[chainNum];
	
    def getX(self, chainNum):
	return self.chain[chainNum].getX();
    def getY(self, chainNum):
	return self.chain[chainNum].getY();
    def getZ(self, chainNum):
	return self.chain[chainNum].getZ();