import matplotlib.pyplot as plt
import Polymer
import EqualAxes
import sys
from matplotlib.animation import FuncAnimation


def update(i):
	plt.cla();
	print('Chain %s has %i atoms.' % (i,polymer.getN(i)));
#	eqAx.push(polymer.getX(i),polymer.getY(i),polymer.getZ(i));
    
    #axMaxRange=eqAx.findMaxRange();
	polymer.plot(i,ax,200);
#	eqAx.push(polymer.getX(i),polymer.getY(i),polymer.getZ(i));
	eqAx.set(ax);

def init():
#    line.set_data([], [])
    eqAx.push(polymer.getX(0),polymer.getY(0),polymer.getZ(0));
    eqAx.set(ax);
    return (None, ax)

if(len(sys.argv)<2):
    print("Error: Give me name of the file you want me to open.");
    print("Hint: You can also pass me number(s) of configuration(s) to plot.\n");
fileName = sys.argv[1];
polymer = Polymer.Polymer(fileName);



fig = plt.figure()
ax = fig.gca(projection='3d');
ax.set_aspect('equal');
eqAx = EqualAxes.EqualAxes();
#eqAx.push(polymer.getX(0),polymer.getY(0),polymer.getZ(0));


	
	
    
#eqAx.set(ax);
anim = FuncAnimation(fig, update, frames=10, interval=200, init_func = init)
plt.show();