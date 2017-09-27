"""
This file should be run from Python shell.
Do it like this:

>>> import Polymer
>>> import sys
>>> polymer = Polymer.Polymer('5dn7.pca')
>>> sys.argv = ['plt.py','50', '80', '110']
>>> execfile('plt.py')

for exit:
>>> exit()

You can close the figure and try other
configurations, by repeating 2 last steps.

you can define all parameters with printung 'set' before the name:
>>> setLineColor = '#000000'
>>> setHueH = 0.5

or return the default parameters:
>>> setLineColor = None
>>> setHueH = 0.05

for exit python:
>>> exit()
"""
#============== default ================
#============ parameters ===============
# Plot
dotSize = None; # if None then optimal size will be found
lineSize = None; # = '#ee0000'; if None then optimal size will be found

dotColor = None; #if None then will be random
lineColor = None; #if None then will be random

dotHueD = 0.05; #[0,1];
dotSaturationD = 0.1; #[0,1];
dotVolumeD = 0.1; #[0,1];

# Axes
elevation = None;
azimut = None;
axisOnOff ='off';
#=======================================

if 'setDotSize' in globals():
    dotSize = setDotSize;
if 'setLineSize' in globals():
    lineSize = setLineSize;

if 'setDotColor' in globals():
    dotColor = setDotColor;
if 'setLineColor' in globals():
    lineColor = setLineColor;

if 'setDotHueD!=None' in globals():
    dotHueD = setDotHueD
if 'setDotSaturationD!=None' in globals():
    dotSaturationD= setDotSaturationD;
if 'setDotVolumeD!=None' in globals():
    dotVolumeD = setDotVolumeD;

import matplotlib.pyplot as plt
import Polymer
import EqualAxes
import Color
import sys

#if(fileName == None):
#    fileName = "1abs_configurations.dat";
#if(polymer==None):
#    polymer = Polymer.Polymer(fileName);

fig = plt.figure()
ax = fig.gca(projection='3d');
ax.set_aspect('equal');
eqAx = EqualAxes.EqualAxes(ax);

if(len(sys.argv)<2):
    confNum = 0;
    eqAx.push(polymer.getX(confNum),polymer.getY(confNum),polymer.getZ(confNum));
    dotSmartColors = Color.arrayWithSmartColors(polymer.getChainLenght(0),
	    dotHueD, dotSaturationD, dotVolumeD, dotColor);
    polymer.plot(confNum, eqAx, dotSize, lineSize, dotSmartColors, lineColor);
    
else:
    for i in range(1,len(sys.argv)):
	confNum = int(sys.argv[i]);
	print('chain %s has %i atoms.' % (sys.argv[i], polymer.getN(confNum)));
	eqAx.push(polymer.getX(confNum),polymer.getY(confNum),polymer.getZ(confNum));
        
    for i in range(1,len(sys.argv)):
	confNum = int(sys.argv[i]);
	dotSmartColors = Color.arrayWithSmartColors(polymer.getChainLenght(0),
		dotHueD, dotSaturationD, dotVolumeD, dotColor);
	polymer.plot(confNum, eqAx, dotSize, lineSize, dotSmartColors, lineColor);
	
    
eqAx.set();
plt.axis(axisOnOff);
plt.show();