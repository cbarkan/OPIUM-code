import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from matplotlib.colors import LogNorm

def cube_show_slider(cube):
    """
    Display a 3d ndarray with a slider to move along the third dimension.

    Extra keyword arguments are passed to imshow
    """
    # check dim
    if not cube.ndim == 3:
        raise ValueError("cube should be an ndarray with ndim == 3")

    # generate figure
    fig = plt.figure()
    ax = plt.subplot(111)
    fig.subplots_adjust(left=0.25, bottom=0.25)

    # select first image
    s = [slice(0, 1) if i == 0 else slice(None) for i in xrange(3)]
    im = cube[s].squeeze()

    # display image
    l = ax.imshow(im,interpolation='none', origin = 'lower',extent=[lower_ext1,upper_ext1,lower_ext0,upper_ext0], norm=LogNorm(vmin=min_value, vmax=max_value))
    cb = fig.colorbar(l,orientation='horizontal')
    plt.ylabel('Coefficient 1')
    plt.xlabel('Coefficient 2')
    # define slider
    axcolor = 'lightgoldenrodyellow'
    ax = fig.add_axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
    
    slider = Slider(ax, 'Coefficient 0', lower_ext2, upper_ext2,
                    valinit=0, valfmt='%5.2f')

    def update(val):
        ind = (int(slider.val) - lower_ext2)/mesh_step
        s = [slice(ind, ind + 1) if i == 0 else slice(None)
                 for i in xrange(3)]
        im = cube[s].squeeze()
        l.set_data(im)

        fig.canvas.draw()

    slider.on_changed(update)

    plt.show()
#_________________________________________________
filename = 'test1'
#WARNING: I think there's an error on line 26 with the bounds of extent. Shouldn't it be ext of coeffs 2 and 1? Also on line 34, it should be lower_ext0, etc, right?

c0_index = np.load(filename + '_c0_index.npy')
c1_index = np.load(filename + '_c1_index.npy')
c2_index = np.load(filename + '_c2_index.npy')

lower_ext0 = c0_index[0]
upper_ext0 = c0_index[-1]
lower_ext1 = c1_index[0]
upper_ext1 = c1_index[-1]
lower_ext2 = c2_index[0]
upper_ext2 = c2_index[-1]
mesh_step = 1.0

mesh = np.load(filename + '_Eig1.npy')

min_value = 1
max_value = 4000

cube_show_slider(mesh)
