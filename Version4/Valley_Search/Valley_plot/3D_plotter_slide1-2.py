import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from matplotlib.colors import LogNorm

def cube_show_slider(cube):
    # check dim
    if not cube.ndim == 3:
        raise ValueError("cube should be an ndarray with ndim == 3")

    # generate figure
    fig = plt.figure()
    ax = plt.subplot(111)
    fig.subplots_adjust(left=0.25, bottom=0.25)

    # select first image
    #s = [slice(0, 1) if i == 0 else slice(None) for i in xrange(3)]
    im = cube[0,:,:]

    # display image
    l = ax.imshow(im,interpolation='none', origin = 'lower',extent=[lower_ext2,upper_ext2,lower_ext1,upper_ext1], norm=LogNorm(vmin=min_value, vmax=max_value))
    cb = fig.colorbar(l,orientation='horizontal')
    plt.ylabel('Coefficient 1')
    plt.xlabel('Coefficient 2')
    # define slider
    axcolor = 'lightgoldenrodyellow'
    ax = fig.add_axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
    
    slider = Slider(ax, 'Coefficient 0', lower_ext0, upper_ext0,
                    valinit=0, valfmt='%5.2f')

    def update(val):
        ind = (slider.val - lower_ext0)/mesh_step
        #s = [slice(ind, ind + 1) if i == 0 else slice(None)
        #         for i in xrange(3)]
        im = cube[ind,:,:]
        l.set_data(im)

        fig.canvas.draw()

    slider.on_changed(update)

    plt.show()
#_________________________________________________
filename = 'startn4-n8-n4-0'

v0_index = np.load(filename + '_v0_index.npy')
v1_index = np.load(filename + '_v1_index.npy')
v2_index = np.load(filename + '_v1_index.npy')
print v0_index
print v1_index
print v2_index

lower_ext0 = v0_index[0]
upper_ext0 = v0_index[-1]
lower_ext1 = v1_index[0]
upper_ext1 = v1_index[-1]
lower_ext2 = v2_index[0]
upper_ext2 = v2_index[-1]
mesh_step = 0.1

mesh = np.load(filename + '_Eig1.npy') + np.load(filename + '_Norm1.npy')

min_value = 4
max_value = 50

cube_show_slider(mesh)
