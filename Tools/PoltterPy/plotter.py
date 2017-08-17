import matplotlib.pyplot as plt
import ChainFile
import EqualAxes
import sys


if(len(sys.argv)<2):
    print("Error: Give me name of the file you want me to open.");
    print("Hint: You can also pass me number(s) of configuration(s) to plot.\n");
fileName = sys.argv[1];
chain = ChainFile.ChainFile(fileName);

fig = plt.figure()
ax = fig.gca(projection='3d');
ax.set_aspect('equal');
eqAx = EqualAxes.EqualAxes();

if(len(sys.argv)<3):
    confNum = 0;
    chain.plot(confNum, ax);
    eqAx.push(chain.getX(confNum),chain.getY(confNum),chain.getZ(confNum));
else:
    for i in range(2,len(sys.argv)):
	confNum = int(sys.argv[i]);
	print('Chain %s has %i atoms.' % (sys.argv[i], chain.getN(confNum)));
	eqAx.push(chain.getX(confNum),chain.getY(confNum),chain.getZ(confNum));
    
    axMaxRange=eqAx.findMaxRange();
    
    for i in range(2,len(sys.argv)):
	confNum = int(sys.argv[i]);
	chain.plot(confNum,ax,axMaxRange);
#	chain.plotOld(confNum, ax);
	
    
eqAx.set(ax);
plt.show();