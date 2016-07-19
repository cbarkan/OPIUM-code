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
    s = [slice(0, 1) if i == 2 else slice(None) for i in xrange(3)]
    im = cube[s].squeeze()

    # display image
    l = ax.imshow(im,interpolation='none', origin = 'lower',extent=[lower_ext1,upper_ext1,lower_ext0,upper_ext0], norm=LogNorm(vmin=min_value, vmax=max_value))
    cb = fig.colorbar(l,orientation='horizontal')
    plt.ylabel('Coefficient 0')
    plt.xlabel('Coefficient 1')
    #ax.annotate('aapoaidscasdlckj',xy=(1,2),xytext = (1,2))
    plt.text(-5,-28,'yuikyuik')
    # define slider
    axcolor = 'lightgoldenrodyellow'
    ax = fig.add_axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
    
    slider = Slider(ax, 'Coefficient 2', lower_ext2, upper_ext2,
                    valinit=0, valfmt='%5.2f')
    #plt.text(0,-5,'yuikyuik')
    def update(val):
        ind = (int(slider.val) - lower_ext2)/mesh_step
        s = [slice(ind, ind + 1) if i == 2 else slice(None)
                 for i in xrange(3)]
        im = cube[s].squeeze()
        l.set_data(im)
        plt.text(-5,-28,'123')

        fig.canvas.draw()

    slider.on_changed(update)

    plt.show()
#_________________________________________________
filename = 'test1'

c0_index = np.load(filename + '_c0_index.npy')
c1_index = np.load(filename + '_c1_index.npy')
c2_index = np.load(filename + '_c2_index.npy')

lower_ext0 = c0_index[0]
upper_ext0 = c0_index[-1]
lower_ext1 = c1_index[0]
upper_ext1 = c1_index[-1]
lower_ext2 = c2_index[0]
upper_ext2 = c2_index[-1]
mesh_step = 5
min_value = 1
max_value = 1000

mesh = np.load(filename + '_Eig1.npy')
cube_show_slider(mesh)
