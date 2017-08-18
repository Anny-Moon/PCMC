"""


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
fileName = sys.argv[1];
polymer = Polymer.Polymer(fileName);

#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

fig = plt.figure()
ax = fig.gca(projection='3d');
ax.set_aspect('equal');
eqAx = EqualAxes.EqualAxes();

#first frame
eqAx.push(polymer.getX(0),polymer.getY(0),polymer.getZ(0));
eqAx.set(ax);
#polymer.plot(0,ax,100);

increment = 5;
#eqAx.set(ax);
anim = animation.FuncAnimation(fig, update,
		    frames=10,
		    interval=200,
		    fargs = (increment,),
		    repeat = False
		    );
#anim.save('movie.gif', dpi=80, writer='imagemagick');
mywriter = animation.FFMpegWriter(fps=60)
anim.save('movie.mp4', writer=mywriter);
plt.show();