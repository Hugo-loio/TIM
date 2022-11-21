import helper as hp
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

plt.style.use('science')

def ldos(fname, name, show: bool):
    data = hp.readfile(fname)
    x = data[0::4,0] + 1
    y = data[1::4,0] + 1
    z = data[2::4,0] + 1
    dos = np.average(data[3::4], axis = 1)
    
    vmax = np.max(dos)
    vmax = vmax if vmax > 0.01 else 1

    fig = plt.figure()
    ax = fig.add_subplot(projection = '3d')
    img = ax.scatter(x, y, z, c = dos, s = 5, vmin = 0, vmax = vmax)

    ax.set(xlabel = r'$x$', ylabel = r'$y$', zlabel = r'$z$')

    '''
    ax.set_xticks(np.arange(0, x[-1] + 1, 5))
    ax.set_yticks(np.arange(0, y[-1] + 1, 5))
    ax.set_zticks(np.arange(0, z[-1] + 1, 5))
    '''

    cbar = fig.colorbar(img, shrink = 0.4)
    cbar.ax.ticklabel_format(useOffset=False)

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 300, bbox_inches='tight')
    if(show):
        plt.show()

def plots(fileNames, names, show):
    for i in range(len(fileNames)):
        ldos(fileNames[i] + ".dat", names[i], show)

#weights = ['2', '3', '4']
weights = ['3']
fileNames = ['ldosBBH3D_L20_w' + w + '_E0_nMu1024_m1.1' for w in weights]
#plots(fileNames, fileNames, False)

#weights = ['2', '3', '4']
weights = ['3']
fileNames = ['ldosBBH3D_L4_w' + w + '_E0_nMu1024_m1.1' for w in weights]
plots(fileNames, fileNames, False)
