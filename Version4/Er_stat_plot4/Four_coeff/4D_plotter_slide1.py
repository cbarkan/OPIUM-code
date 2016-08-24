import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from matplotlib.colors import LogNorm

def cube_show_slider(cube):
    # check dim
    if not cube.ndim == 4:
        raise ValueError("cube should be an ndarray with ndim == 4")

    # generate figure
    fig = plt.figure()
    ax = plt.subplot(111)
    fig.subplots_adjust(left=0.25, bottom=0.25)

    # select first image
    #s = [slice(0, 0, 1) if i == 0 else slice(None) for i in xrange(4)]
    #im = cube[0,0,:,:].squeeze()
    im = cube[0,0,:,:]

    # display image
    l = ax.imshow(im,interpolation='none', origin = 'lower',extent=[lower_ext3,upper_ext3,lower_ext2,upper_ext2], norm=LogNorm(vmin=min_value, vmax=max_value))
    cb = fig.colorbar(l,orientation='horizontal')
    plt.ylabel('Coefficient 2')
    plt.xlabel('Coefficient 3')
    # define slider
    axcolor = 'lightgoldenrodyellow'
    coeff0_ax = fig.add_axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
    coeff1_ax = fig.add_axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)
    
    slider0 = Slider(coeff0_ax, 'Coefficient 0', lower_ext0, upper_ext0,
                    valinit=lower_ext0, valfmt='%5.2f')
    
    slider1 = Slider(coeff1_ax, 'Coefficient 1', lower_ext1, upper_ext1,
                    valinit=lower_ext1, valfmt='%5.2f')

    def update(val):
        ind0 = (slider0.val - lower_ext0)/mesh_step
        ind1 = (slider1.val - lower_ext1)/mesh_step
        #s = [slice(ind, ind + 1) if i == 0 else slice(None)
        #         for i in xrange(4)]
        im = cube[ind0,ind1,:,:]
        l.set_data(im)

        fig.canvas.draw()

    slider0.on_changed(update)
    slider1.on_changed(update)

    plt.show()
#_________________________________________________
filename = 'vs6_result'

c0_index = np.load(filename + '_c0_index.npy')
c1_index = np.load(filename + '_c1_index.npy')
c2_index = np.load(filename + '_c2_index.npy')
c3_index = np.load(filename + '_c3_index.npy')

lower_ext0 = c0_index[0]
upper_ext0 = c0_index[-1]
lower_ext1 = c1_index[0]
upper_ext1 = c1_index[-1]
lower_ext2 = c2_index[0]
upper_ext2 = c2_index[-1]
lower_ext3 = c3_index[0]
upper_ext3 = c3_index[-1]
mesh_step = 0.0000005/7.0

print 'lower_ext0 = %s' % lower_ext0
print 'lower_ext1 = %s' % lower_ext1

mesh = np.load(filename + '_Eig1.npy') + np.load(filename + '_Norm1.npy')
#mesh = mesh_real*10000.0
min_value = np.min(mesh)
max_value = np.max(mesh)
print min_value
print max_value

cube_show_slider(mesh)
