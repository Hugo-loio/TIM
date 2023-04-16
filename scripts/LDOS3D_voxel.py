import helper as hp
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import gridspec
import numpy as np

#plt.style.use('science')
plt.rcParams.update({
    "text.usetex": True,
    #"font.family": "Helvetica"
    "font.size" : 16
})

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
    #dosColors[:,-1] = np.array(s)
    #dosColors[:,-1] = 0.8
    norm = Normalize(vmin = 0, vmax = vmax)

    spec = gridspec.GridSpec(ncols=2, nrows=1, wspace = 0, width_ratios=[40, 1])
    fig = plt.figure()
    ax = fig.add_subplot(spec[0], projection = '3d')
    #ax.view_init(30, 210)

    voxels = np.ones((int(np.amax(x)),int(np.amax(y)),int(np.amax(z))))
    facecolors = [[[(0,0,0,1) for z in y] for y in x] for x in voxels]
    edgecolors = [[[(0,0,0,1) for z in y] for y in x] for x in voxels]
    for i,c in enumerate(dosColors):
        facecolors[int(x[i]-1)][int(y[i]-1)][int(z[i]-1)] = dosColors[i]
        edgecolors[int(x[i]-1)][int(y[i]-1)][int(z[i]-1)] = (0,0,0,dosColors[i,-1])
    ax.voxels(voxels, edgecolors = np.array(edgecolors), facecolors = np.array(facecolors), shade = False)

    #fig.set_size_inches(4, 3)
    ax.set(xlabel = r'$x$', ylabel = r'$y$', zlabel = r'$z$')
    ax.set_xlabel(r'$x$', labelpad = -4, fontsize = 20)
    ax.set_ylabel(r'$y$', labelpad = -4, fontsize = 20)
    #ax.tick_params(axis='y', which='minor', pad=0)
    ax.set_zlabel(r'$z$', labelpad = -8, fontsize = 20)
    ax.invert_xaxis()
    ax.invert_zaxis()

    xticks = ax.get_xticks()[ax.get_xticks().tolist().index(0):-2] + 0.5
    xlabels = [str(int(tick + 0.5)) for tick in xticks]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)
    yticks = ax.get_yticks()[ax.get_yticks().tolist().index(0):-2] + 0.5
    ylabels = [str(int(tick + 0.5)) for tick in yticks]
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)
    zticks = ax.get_zticks()[ax.get_zticks().tolist().index(0):-1] + 0.5
    zlabels = [str(int(tick + 0.5)) for tick in zticks]
    ax.set_zticks(zticks)
    ax.set_zticklabels(zlabels)

    ax.tick_params(axis='both', which='major', pad=-3)
    ax2 = fig.add_subplot(spec[1])
    ax2.set_aspect(20/vmax)
    pos2 = ax2.get_position()
    ax2.set_position([0.95*pos2.x0, pos2.y0, pos2.width, pos2.height])
#pos2 = [pos1.x0 + 0.3, pos1.y0 + 0.3,  pos1.width / 2.0, pos1.height / 2.0]
    cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax = ax2, orientation = 'vertical')
    cbar.ax.set_title(r'$\rho(0, \textbf{r})$')

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 300, bbox_inches='tight', pad_inches = 0)
    fig.savefig(hp.plot_dir() + name + ".pdf", bbox_inches = 'tight', pad_inches = 0)
    if(show):
        plt.show()

def plots(fileNames, names, show):
    for i in range(len(fileNames)):
        ldos(fileNames[i] + ".dat", names[i], show)

weights = ['3', '4']
#sizes = ['10', '16', '20', '24', '30']
sizes = ['30']
#sizes = sizes[0:1]
fileNames = ['ldosBBH3D_L' + size + '_w' + w + '_E0_nMu1024_m1.1' for w in weights for size in sizes]
names = [name + '_voxel' for name in fileNames]
plots(fileNames, names, False)
