programDescription ="""\
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
fps = 2;
dpi = 200;
frames = None; #number of frames, defailt = all frames

# Plot
dotSize = 8;
dotColor = '#cc0000';
lineColor = '#006699';

# Axes
elevation = None;
azimut = None;
axisOnOff ='off';
#=======================================

import matplotlib.pyplot as plt
import Polymer
import EqualAxes
import sys
import math
import matplotlib.animation as animation


def update(i, increment):
	confNum = increment+(i-1)*increment;
	plt.cla();
	print('Chain %s has %i atoms.' % (confNum,polymer.getN(confNum)));
    
    #axMaxRange=eqAx.findMaxRange();
	polymer.plot(confNum,ax,800/dotSize, dotColor, lineColor);
	eqAx.set(ax);
	ax.view_init(elevation, azimut);
	plt.axis(axisOnOff);

if(len(sys.argv)<2):
    print(programDescription);
    exit();

fileNameIn = sys.argv[1];
polymer = Polymer.Polymer(fileNameIn);

fig = plt.figure()
ax = fig.gca(projection='3d');
ax.set_aspect('equal');
eqAx = EqualAxes.EqualAxes();

#first frame
eqAx.push(polymer.getX(0),polymer.getY(0),polymer.getZ(0));
eqAx.set(ax);
ax.view_init(elevation, azimut);
plt.axis(axisOnOff);

increment = 1;

#if pass increment
if(len(sys.argv)>2):
    increment = int(sys.argv[2]);
    if(increment < 1):
	increment = 1;
if(frames==None):
    frames =int(math.ceil(polymer.getNumChains()/float(increment)));
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
    fileNameOut=fileName + extention;
    
    #if pass fileNameOut
    if(len(sys.argv)>4):
	fileNameOut = sys.argv[4] + extention;
    
    if(sys.argv[3]=='gif'):
	writer = animation.ImageMagickFileWriter(fps=fps);
    if(sys.argv[3]=='mp4'):
	writer = animation.FFMpegWriter(fps=fps);
	
    anim.save(fileNameOut, writer=writer, dpi=dpi);
    print("File %s was successfully generated!" % fileNameOut);
