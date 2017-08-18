"""
python ./movieMaker.py <path/fileIn.ext> <increment> <gif/mp4> <path/fileOut>
python ./movieMaker.py 5dn7.pca 3 gif movie

For saving option you need additional packages which
are not included in 'matplotlib'.

For saving gif 'imagemagick' pack is needed.
For saving mp4 'ffmpeg' pack is needed.
"""
import matplotlib.pyplot as plt
import Polymer
import EqualAxes
import sys
import matplotlib.animation as animation

def update(i, increment):
	confNum = increment+(i-1)*increment;
	plt.cla();
	print('Chain %s has %i atoms.' % (confNum,polymer.getN(confNum)));
    
    #axMaxRange=eqAx.findMaxRange();
	polymer.plot(confNum,ax,None,'#e60000', '#006699');
	eqAx.set(ax);

def init():
#    line.set_data([], [])
    eqAx.push(polymer.getX(0),polymer.getY(0),polymer.getZ(0));
    eqAx.set(ax);
    polymer.plot(0,ax,None,'#000000', '#000000');
    return (ax)

if(len(sys.argv)<2):
    print("Error: Give me name of the file you want me to open.");
    print("Hint: You can also pass me number(s) of configuration(s) to plot.\n");
fileNameIn = sys.argv[1];
polymer = Polymer.Polymer(fileNameIn);

fig = plt.figure()
ax = fig.gca(projection='3d');
ax.set_aspect('equal');
eqAx = EqualAxes.EqualAxes();

#first frame
eqAx.push(polymer.getX(0),polymer.getY(0),polymer.getZ(0));
eqAx.set(ax);

increment = 1;

#if pass increment
if(len(sys.argv)>2):
    increment = int(sys.argv[2]);
    if(increment < 1):
	increment = 1;

anim = animation.FuncAnimation(fig, update,
		    frames=3,
		    interval=200,
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
	anim.save('movie.gif', dpi=80, writer='imagemagick');
    if(sys.argv[3]=='mp4'):
	writer = animation.FFMpegWriter(fps=2);
	anim.save('movie.mp4', writer=writer);
    print("File %s was successfully generated!" % fileNameOut);
