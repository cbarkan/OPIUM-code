from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider, Button, RadioButtons

"""
def update(val):
    v0 = slider0.val
    l.set_ydata(amp*np.sin(2*np.pi*freq*t))
    fig.canvas.draw_idle()
sfreq.on_changed(update)
samp.on_changed(update)
"""

filename = 'startn4-n8-n4-0'

v0_index = np.load(filename + '_v0_index.npy')
v1_index = np.load(filename + '_v1_index.npy')
v2_index = np.load(filename + '_v1_index.npy')

lower_ext0 = v0_index[0]
upper_ext0 = v0_index[-1]
lower_ext1 = v1_index[0]
upper_ext1 = v1_index[-1]
lower_ext2 = v2_index[0]
upper_ext2 = v2_index[-1]
mesh_step = 0.1

mesh = np.load(filename + '_Eig1.npy') + np.load(filename + '_Norm1.npy')

min_value = 1
max_value = 50

fig = plt.figure()
ax = fig.gca(projection='3d')
X = v1_index.copy()
Y = v2_index.copy()
X, Y = np.meshgrid(X, Y)
R = mesh[1,:,:]
surf = ax.plot_surface(X, Y, R, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
ax.set_zlim(min_value, max_value)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=0.5, aspect=5)

#slider0 = Slider(v0_val)

plt.show()



