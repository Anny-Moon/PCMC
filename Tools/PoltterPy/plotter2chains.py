import matplotlib.pyplot as plt
import Polymer
import EqualAxes
import sys

fileName = "0confR.dat";
polymer = Polymer.Polymer(fileName);

fileName2 = "0confR2.dat";
polymer2 = Polymer.Polymer(fileName2);

fig = plt.figure()
ax = fig.gca(projection='3d');
ax.set_aspect('equal');
eqAx = EqualAxes.EqualAxes();

if(len(sys.argv)<2):
    confNum = 0;
    polymer.plot(confNum,ax,None,'#e60000','#006699');
    polymer2.plot(confNum,ax,None,'#006699','#e60000');
    eqAx.push(polymer.getX(confNum),polymer.getY(confNum),polymer.getZ(confNum));
    eqAx.push(polymer2.getX(confNum),polymer2.getY(confNum),polymer2.getZ(confNum));
else:
    for i in range(1,len(sys.argv)):
	confNum = int(sys.argv[i]);
	print('Chain %s has %i atoms.' % (sys.argv[i], polymer.getN(confNum)));
	eqAx.push(polymer.getX(confNum),polymer.getY(confNum),polymer.getZ(confNum));
	eqAx.push(polymer2.getX(confNum),polymer2.getY(confNum),polymer2.getZ(confNum));
    
    axMaxRange=eqAx.findMaxRange();
    
    for i in range(1,len(sys.argv)):
	confNum = int(sys.argv[i]);
	polymer.plot(confNum,ax,axMaxRange,'#e60000','#006699');
	polymer2.plot(confNum,ax,axMaxRange,'#006699','#e60000');
#	polymer.plotOld(confNum, ax);
	
    
eqAx.set(ax);
plt.show();