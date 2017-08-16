import matplotlib.pyplot as plt
import ChainFile
import EqualAxes
import sys

fileName = "0confR.dat";
chain = ChainFile.ChainFile(fileName);

fileName2 = "0confR2.dat";
chain2 = ChainFile.ChainFile(fileName2);

fig = plt.figure()
ax = fig.gca(projection='3d');
ax.set_aspect('equal');
eqAx = EqualAxes.EqualAxes();

if(len(sys.argv)<2):
    confNum = 0;
    chain.plot(confNum,ax,None,'#e60000','#006699');
    chain2.plot(confNum,ax,None,'#006699','#e60000');
    eqAx.push(chain.getX(confNum),chain.getY(confNum),chain.getZ(confNum));
    eqAx.push(chain2.getX(confNum),chain2.getY(confNum),chain2.getZ(confNum));
else:
    confNum = int(sys.argv[1]);
    print('Chain %s has %i atoms.' % (sys.argv[1], chain.getN(confNum)));
    eqAx.push(chain.getX(confNum),chain.getY(confNum),chain.getZ(confNum));
    eqAx.push(chain2.getX(confNum),chain2.getY(confNum),chain2.getZ(confNum));
    
    axMaxRange=eqAx.findMaxRange();
    
    for i in range(1,len(sys.argv)):
	confNum = int(sys.argv[i]);
	chain.plot(confNum,ax,axMaxRange,'#e60000','#006699');
	chain2.plot(confNum,ax,axMaxRange,'#006699','#e60000');
#	chain.plotOld(confNum, ax);
	
    
eqAx.set(ax);
plt.show();