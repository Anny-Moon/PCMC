"""
Run me like this:
python ./plotter.py <path/fileIn.ext> <configurations>
python ./movieMaker.py data/5dn7_rescaled.pca 30 44 58

If you pass only the first argument then I will
plot the  0-th configuration.
"""

#============ parameters ===============
# Plot
randomSeed = 11;

dotSize = 3.3; #if None then optimal size will be found
lineSize = 0.5; # = '#ee0000'; if None then optimal size will be found

#dotColor = (0.0, 0.5, 0.8); #if None then will be random
#lineColor = (0.0, 0.4, 0.7); #if None then will be random
#dotColor2 = (0.8, 0.0, 0.2); #if None then will be random
#lineColor2 = (0.7, 0.0, 0.1); #if None then will be random

#lineColor = (17./255.,33./255.,130./255.);
#dotColor = (47./255.,116./255.,161./255.);

dotColor = (244, 0, 65);
dotColor = tuple(float(x)/255. for x in dotColor);

lineColor = (137, 0, 36);
lineColor = tuple(float(x)/255. for x in lineColor);

dotColor2 = (50, 135, 210);
dotColor2 = tuple(float(x)/255. for x in dotColor2);

lineColor2 = (15, 61, 101);
lineColor2 = tuple(float(x)/255. for x in lineColor2);

backgroundColor = (255, 255, 255);
backgroundColor = tuple(float(x)/255. for x in backgroundColor);

#backgroundColor = (255./255.,255./255.,255./255.);

# Axes
elevation = 0;
azimut = 300;
axisOnOff ='on';
#=======================================


import MayaPolymer
import sys
from mayavi import mlab


if(len(sys.argv)<2):
    print(__doc__);
    exit();
#fileName = sys.argv[1];
fileName = "/Users/annsi118/Documents/Git_projects/PCMC/Projects/MonteCarlo2chains/results/Configurations/0confR.dat";
fileName2 = "/Users/annsi118/Documents/Git_projects/PCMC/Projects/MonteCarlo2chains/results/Configurations/0confR2.dat";

#fileName = "0confR.dat";
#fileName2 = "0confR2.dat";

polymer = MayaPolymer.Polymer(fileName);
polymer2 = MayaPolymer.Polymer(fileName2);

fig=mlab.figure(bgcolor = backgroundColor,size=(1000,1000))

if(len(sys.argv)<3):
    confNum = 0;
    polymer.plot(confNum, fig, dotSize, lineSize, dotColor, lineColor);
    polymer2.plot(confNum, fig, dotSize, lineSize, dotColor2, lineColor2);
    
else:

    for i in range(2,len(sys.argv)):
	confNum = int(sys.argv[i]);
	print('Chain %s has %i atoms.' % (sys.argv[i], polymer.getN(confNum)));
	
	polymer.plot(confNum, fig, dotSize, lineSize, dotColor, lineColor);
	polymer2.plot(confNum, fig, dotSize, lineSize, dotColor2, lineColor2);


# Your fantastic plotting script

#mlab.save('1.png')

mlab.show()

