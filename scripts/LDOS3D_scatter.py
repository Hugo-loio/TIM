import helper as hp
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import gridspec
import numpy as np

#plt.style.use('science')

def ldos(fname, name, show: bool):
    data = hp.readfile(fname)
    x = data[0::4,0] + 1
    y = data[1::4,0] + 1
    z = data[2::4,0] + 1
    #dos = np.average(data[3::4], axis = 1) - 0.4
    dos = np.average(data[3::4], axis = 1)

    s = []
    maxD = np.sqrt(3*np.amax(x)**2)
    for i in range(len(x)):
        s.append(np.exp(-3*(np.sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]) - np.sqrt(3))/maxD))
    
    vmax = np.max(dos)
    vmax = vmax if vmax > 0.01 else 1

    colors = ["palegreen", "maroon"]
    nodes = [0.0, 1.0]
    cmap = LinearSegmentedColormap.from_list("mycmap", list(zip(nodes,colors)))
    dosColors = cmap(dos/vmax)
    dosColors[:,-1] = np.array(s)
    norm = Normalize(vmin = 0, vmax = vmax)

    spec = gridspec.GridSpec(ncols=2, nrows=1, width_ratios=[40, 1])
    fig = plt.figure()
    ax = fig.add_subplot(spec[0], projection = '3d')
    #ax.view_init(30, 210)
    img = ax.scatter(x, y, z, c = dosColors, s = 100)

    ax.set(xlabel = r'$x$', ylabel = r'$y$', zlabel = r'$z$')
    ax.set_zlabel(r'$z$', labelpad = 0)
    ax.invert_xaxis()
    ax.invert_zaxis()

    '''
    ax.set_xticks(np.arange(0, x[-1] + 1, 5))
    ax.set_yticks(np.arange(0, y[-1] + 1, 5))
    ax.set_zticks(np.arange(0, z[-1] + 1, 5))
    '''
    ax.tick_params(axis='both', which='major', pad=0)
    ax2 = fig.add_subplot(spec[1])

    #divider = make_axes_locatable(ax)
    #cax = divider.append_axes('right', size='5%', pad=0.05)

    ax2.set_aspect(20/vmax)
    cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax = ax2, orientation = 'vertical')
    cbar.ax.set_title(r'$\rho(0, \overrightarrow{n})$')
    #cbar.ax.ticklabel_format(useOffset=False)

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 300, bbox_inches='tight', pad_inches = 0)
    if(show):
        plt.show()

def plots(fileNames, names, show):
    for i in range(len(fileNames)):
        ldos(fileNames[i] + ".dat", names[i], show)

weights = ['2', '3', '4']
fileNames = ['ldosBBH3D_L20_w' + w + '_E0_nMu1024_m1.1' for w in weights]
#plots(fileNames, fileNames, False)

weights = ['2', '3', '4']
fileNames = ['ldosBBH3D_L20_w' + w + '_E0_nMu1024_m1.1' for w in weights]
#plots(fileNames, fileNames, False)

weights = ['2', '3', '4']
fileNames = ['ldosBBH3D_L10_w' + w + '_E0_nMu1024_m1.1' for w in weights]
#plots(fileNames, fileNames, False)

weights = ['3', '4']
sizes = ['10', '16', '20', '24']
fileNames = ['ldosBBH3D_L' + size + '_w' + w + '_E0_nMu1024_m1.1' for w in weights for size in sizes]
plots(fileNames, fileNames, False)


#weights = ['2', '3', '4']
weights = ['3']
fileNames = ['ldosBBH3D_L4_w' + w + '_E0_nMu1024_m1.1' for w in weights]
#plots(fileNames, fileNames, False)
