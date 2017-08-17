"""
This file should be run from Python shell.
Do it like this:

>>> import Polymer
>>> import sys
>>> polymer = Polymer.Polymer('5dn7.pca')
>>> sys.argv = ['plt.py','50', '80', '110']
>>> execfile('plt.py')

You can close the figure and try other
configurations, by repeating 2 last steps.
"""

import matplotlib.pyplot as plt
import Polymer
import EqualAxes
import sys

#if(fileName == None):
#    fileName = "1abs_configurations.dat";
#if(polymer==None):
#    polymer = Polymer.Polymer(fileName);

fig = plt.figure()
ax = fig.gca(projection='3d');
ax.set_aspect('equal');
eqAx = EqualAxes.EqualAxes();

if(len(sys.argv)<2):
    confNum = 0;
    polymer.plot(confNum, ax);
    eqAx.push(polymer.getX(confNum),polymer.getY(confNum),polymer.getZ(confNum));
else:
    for i in range(1,len(sys.argv)):
	confNum = int(sys.argv[i]);
	print('chain %s has %i atoms.' % (sys.argv[i], polymer.getN(confNum)));
	eqAx.push(polymer.getX(confNum),polymer.getY(confNum),polymer.getZ(confNum));
    
    axMaxRange=eqAx.findMaxRange();
    
    for i in range(1,len(sys.argv)):
	confNum = int(sys.argv[i]);
	polymer.plot(confNum,ax,axMaxRange);
#	polymer.plotOld(confNum, ax);
	
    
eqAx.set(ax);
plt.show();