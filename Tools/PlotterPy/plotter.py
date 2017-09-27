"""
Run me like this:
python ./plotter.py <path/fileIn.ext> <configurations>
python ./movieMaker.py data/5dn7_rescaled.pca 30 44 58

If you pass only the first argument then I will
plot the  0-th configuration.
"""

#============ parameters ===============
# Plot
dotSize = None; #if None then optimal size will be found
lineSize = None; # = '#ee0000'; if None then optimal size will be found

dotColor = None; #if None then will be random
lineColor = None; #if None then will be random

dotHueDispersion = 0.05; #[0,1];
dotSaturationDispersion = 0.1; #[0,1];
dotVolumeDispersion = 0.1; #[0,1];

# Axes
elevation = None;
azimut = None;
axisOnOff ='off';
#=======================================


import matplotlib.pyplot as plt
import Polymer
import EqualAxes
import Color
import sys


if(len(sys.argv)<2):
    print(__doc__);
    exit();
fileName = sys.argv[1];
polymer = Polymer.Polymer(fileName);

fig = plt.figure()
ax = fig.gca(projection='3d');
ax.set_aspect('equal');
eqAx = EqualAxes.EqualAxes(ax);

if(len(sys.argv)<3):
    confNum = 0;
    eqAx.push(polymer.getX(confNum),polymer.getY(confNum),polymer.getZ(confNum));
    dotSmartColors = Color.arrayWithSmartColors(polymer.getChainLenght(0),
	    dotHueDispersion, dotSaturationDispersion, dotVolumeDispersion, dotColor);
    polymer.plot(confNum, eqAx, dotSize, lineSize, dotSmartColors, lineColor);
    
else:
    for i in range(2,len(sys.argv)):
	confNum = int(sys.argv[i]);
	print('Chain %s has %i atoms.' % (sys.argv[i], polymer.getN(confNum)));
	eqAx.push(polymer.getX(confNum),polymer.getY(confNum),polymer.getZ(confNum));

	
    for i in range(2,len(sys.argv)):
	confNum = int(sys.argv[i]);
	dotSmartColors = Color.arrayWithSmartColors(polymer.getChainLenght(0),
		dotHueDispersion, dotSaturationDispersion, dotVolumeDispersion, dotColor);
	polymer.plot(confNum, eqAx, dotSize, lineSize, dotSmartColors, lineColor);
#	polymer.smartColorPlot(confNum,ax, axMaxRange, "#002255");
#	polymer.happyPlot(confNum,ax, axMaxRange);
#	polymer.plotOld(confNum, ax);
	
    
eqAx.set();
plt.axis(axisOnOff);
plt.show();