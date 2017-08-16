import matplotlib.pyplot as plt
import ChainFile
import EqualAxes
import sys

fileName = "0confR.dat";
chain = ChainFile.ChainFile(fileName);

fig = plt.figure()
ax = fig.gca(projection='3d');
ax.set_aspect('equal');
eqAx = EqualAxes.EqualAxes();

if(len(sys.argv)<2):
    confNum = 0;
    chain.plot(confNum, ax, 40);
#    eqAx.push(chain.getChain(confNum).getX(),chain.getChain(confNum).getY(),chain.getChain(confNum).getZ());
    eqAx.push(chain.getX(confNum),chain.getY(confNum),chain.getZ(confNum));
else:
    for i in range(1,len(sys.argv)):
	confNum = int(sys.argv[i]);
	print('Chain %s has %i atoms.' % (sys.argv[i], chain.getN(confNum)));
#	eqAx.push(chain.getChain(confNum).getX(),chain.getChain(confNum).getY(),chain.getChain(confNum).getZ());
	eqAx.push(chain.getX(confNum),chain.getY(confNum),chain.getZ(confNum));
    
    axMaxRange=eqAx.findMaxRange();
    
    for i in range(1,len(sys.argv)):
	confNum = int(sys.argv[i]);
	
	chain.plot(confNum,ax,axMaxRange);
#	chain.plotOld(confNum, ax);
	
    
eqAx.set(ax);
plt.show();