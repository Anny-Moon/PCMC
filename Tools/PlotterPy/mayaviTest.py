import numpy as np
from mayavi import mlab

n_mer, n_long = 6, 11
pi = np.pi
dphi = pi/1000.0
phi = np.arange(0.0, 2*pi + 0.5*dphi, dphi)
mu = phi*n_mer
x = np.cos(mu)*(1+np.cos(n_long*mu/n_mer)*0.5)
y = np.sin(mu)*(1+np.cos(n_long*mu/n_mer)*0.5)
z = np.sin(n_long*mu/n_mer)*0.5

l = mlab.plot3d(x, y, z, np.sin(mu),tube_radius=0.025, colormap='Spectral')
mlab.view(distance=4.75);
mlab.pitch(-2.0)

mlab.show()