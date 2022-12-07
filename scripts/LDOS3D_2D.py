import helper as hp
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

plt.style.use('science')

def ldos(ax, fname, label):
    data = hp.readfile(fname)
    x = data[0::4,0] + 1
    y = data[1::4,0] + 1
    z = data[2::4,0] + 1
    #dos = np.average(data[3::4], axis = 1) - 0.4
    dos = np.average(data[3::4], axis = 1)
    allDist = np.sqrt(x*x + y*y + z*z) - np.sqrt(3)
    d = []
    dDos = []
    for dist in allDist:
        if(not np.any(d == dist)):
            d.append(dist)
            dDos.append(np.average(dos[np.where(allDist == dist)[0]]))
    
    sort = np.argsort(d)
    d = np.array(d)[sort]
    dDos = np.array(dDos)[sort]

    ax.plot(d, dDos, label = label)


def plots(fileNames, labels, show):
    fig, ax = plt.subplots()
    ax.set(xlabel = r'$| \overrightarrow{n} |$', ylabel = r'$\rho(0,| \overrightarrow{n} |)$')

    for i in range(len(fileNames)):
        ldos(ax, fileNames[i] + ".dat", labels[i])

    plt.yscale('log')
    ax.margins(x = 0)
    ax.legend(fontsize = 6, ncol = 1)

    fig.savefig(hp.plot_dir() + "ldosBBH3D_2D_E0_nMu1024_m1.1.png", dpi = 300, bbox_inches='tight')
    if(show):
        plt.show()

weights = ['3', '4']
sizes = ['10', '16', '20','24']
fileNames = ['ldosBBH3D_L' + size + '_w' + w + '_E0_nMu1024_m1.1' for w in weights for size in sizes]
labels = ['W = ' + w + ', L = ' + size for w in weights for size in sizes]
plots(fileNames, labels, False)
