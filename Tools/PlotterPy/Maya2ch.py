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

dotSize = None; #if None then optimal size will be found
lineSize = None; # = '#ee0000'; if None then optimal size will be found

dotColor = (0.0, 0.5, 0.8); #if None then will be random
lineColor = (0.0, 0.2, 0.4); #if None then will be random
dotColor2 = (0.8, 0.0, 0.2); #if None then will be random
lineColor2 = (0.3, 0.0, 0.1); #if None then will be random

dotHueDispersion = 0.05; #[0,1];
dotSaturationDispersion = 0.1; #[0,1];
dotVolumeDispersion = 0.1; #[0,1];

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
polymer = MayaPolymer.Polymer(fileName);
polymer2 = MayaPolymer.Polymer(fileName2);

if(len(sys.argv)<3):
    confNum = 0;
    polymer.plot(confNum, dotSize, lineSize, dotColor, lineColor);
    polymer2.plot(confNum, dotSize, lineSize, dotColor2, lineColor2);
    
else:

    for i in range(2,len(sys.argv)):
	confNum = int(sys.argv[i]);
	print('Chain %s has %i atoms.' % (sys.argv[i], polymer.getN(confNum)));
	
	polymer.plot(confNum, dotSize, lineSize, dotColor, lineColor);
	polymer2.plot(confNum, dotSize, lineSize, dotColor2, lineColor2);


# Your fantastic plotting script

#mlab.save('pure_beauty.png')	
f=mlab.figure(bgcolor = (1, 1, 1))
mlab.show()

