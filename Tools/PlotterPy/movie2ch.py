"""\
Run it like this:
python ./movieMaker.py <path/fileIn.ext> <increment> <gif/mp4> <path/fileOut>
python ./movieMaker.py 5dn7.pca 3 gif movie

If you pass only the first (and second) argument then
no file will be generated, but everything will be 
showen on screen.

!!!
For saving option you need additional packages which
are not included in matplotlib.

For saving gif imagemagick pack is needed.
For saving mp4 ffmpeg pack is needed.
"""
#============ parameters ===============
# Video
fps = 3;
dpi = 250;
frames = None; #number of frames, defailt = all frames

# Plot
dotSize = 20; #if None then optimal size will be found
lineSize = 1; #if None then optimal size will be found

dotColor = "#cc0000"; #if None then will be random
lineColor = "#660000"; #if None then will be random
dotColor2 = "#0066cc"; #if None then will be random
lineColor2 = "#00264d"; #if None then will be random

dotHueDispersion = 0.02; #[0,1];
dotSaturationDispersion = 0.5; #[0,1];
dotVolumeDispersion = 0.2; #[0,1];

# Axes
elevation = None;
azimut = None;
axisOnOff ='on';
#=======================================

import matplotlib.pyplot as plt
import Polymer
import EqualAxes
import Color
import sys
import math
import matplotlib.animation as animation


def update(i, increment):
	confNum = increment+(i-1)*increment;
	plt.cla();
	print('Chain %s has %i atoms.' % (confNum,polymer.getN(confNum)));
    
    #axMaxRange=eqAx.findMaxRange();
	eqAx.clean();
	eqAx.push(polymer.getX(confNum),polymer.getY(confNum),polymer.getZ(confNum));
	eqAx.push(polymer2.getX(confNum),polymer2.getY(confNum),polymer2.getZ(confNum));
	polymer.plot(confNum, eqAx, dotSize, lineSize, dotSmartColors, lineColor);
	polymer2.plot(confNum, eqAx, dotSize, lineSize, dotSmartColors2, lineColor2);
#	polymer.smartColorPlot(confNum,ax,800/dotSize, dotColor, lineColor);
	eqAx.set();
	ax.view_init(elevation, azimut);
	plt.axis(axisOnOff);

if(len(sys.argv)<2):
    print(__doc__);
    exit();

fileNameIn = "test.txt";
fileName = "/Users/annsi118/Documents/Git_projects/PCMC/Projects/MonteCarlo2chains/results/Configurations/0confR.dat";
fileName2 = "/Users/annsi118/Documents/Git_projects/PCMC/Projects/MonteCarlo2chains/results/Configurations/0confR2.dat";

polymer = Polymer.Polymer(fileName);
polymer2 = Polymer.Polymer(fileName2);

dotSmartColors = Color.arrayWithSmartColors(polymer.getChainLenght(0),
		dotHueDispersion, dotSaturationDispersion, dotVolumeDispersion, dotColor);
dotSmartColors2 = Color.arrayWithSmartColors(polymer2.getChainLenght(0),
		dotHueDispersion, dotSaturationDispersion, dotVolumeDispersion, dotColor2);
fig = plt.figure()
ax = fig.gca(projection='3d');
ax.set_aspect('equal');
eqAx = EqualAxes.EqualAxes(ax);

#first frame
eqAx.push(polymer.getX(0),polymer.getY(0),polymer.getZ(0));
eqAx.push(polymer2.getX(0),polymer2.getY(0),polymer2.getZ(0));
eqAx.set();
ax.view_init(elevation, azimut);
plt.axis(axisOnOff);

increment = 1;

#if pass increment
if(len(sys.argv)>2):
    increment = int(sys.argv[2]);
    if(increment < 1):
	increment = 1;
if(frames==None):
    frames =int(math.ceil(polymer2.getNumChains()/float(increment)));
    print("Number of frames: %s."%frames);
anim = animation.FuncAnimation(fig, update,
		    frames=frames,
		    interval=1000/fps,
		    fargs = (increment,),
		    repeat = False
		    );

#if no saving
if(len(sys.argv)<4):
    plt.show();

#if save
else:
    extention = sys.argv[3];
    fileName=fileNameIn[:-3];
    fileNameOut=fileName +extention;
    
    #if pass fileNameOut
    if(len(sys.argv)>4):
	fileNameOut = sys.argv[4] +"."+ extention;
    
    if(sys.argv[3]=='gif'):
	writer = animation.ImageMagickFileWriter(fps=fps);
    if(sys.argv[3]=='mp4'):
	writer = animation.FFMpegWriter(fps=fps);
	
    anim.save(fileNameOut, writer=writer, dpi=dpi);
    print("File %s was successfully generated!" % fileNameOut);
