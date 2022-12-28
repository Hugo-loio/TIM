import helper as hp
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

plt.style.use('science')

def ldos(ax, fname, label):
    data = hp.readfile(fname)
    d = data[0::2,0]
    dos = np.average(data[1::2], axis = 1)
    dos_err = np.std(data[1::2], axis = 1)/np.sqrt(len(data[0]))

    #ax.plot(d, dos, label = label)
    ax.errorbar(d, dos, yerr = dos_err, capsize = 2, linewidth = 0.5, capthick = 0.5, label = label)


def plots(name, fileNames, labels, show):
    fig, ax = plt.subplots()
    ax.set(xlabel = r'$| \overrightarrow{n} |$', ylabel = r'$\rho(0,| \overrightarrow{n} |)$')

    for i in range(len(fileNames)):
        ldos(ax, fileNames[i] + ".dat", labels[i])

    plt.yscale('log')
    #plt.xscale('log')
    #ax.margins(x = 0)
    ax.legend(fontsize = 6, ncol = 1)

    fig.savefig(hp.plot_dir() + name, dpi = 300, bbox_inches='tight')
    if(show):
        plt.show()

weights = ['3', '4']
sizes = ['20', '30', '40']
fileNames = ['ldosBBH3D_L' + size + '_w' + w + '_E0_nMu1024_m1.1_edge' for w in weights for size in sizes]
labels = ['W = ' + w + ', L = ' + size for w in weights for size in sizes]
plots("ldosBBH3D_2D_E0_nMu1024_edge", fileNames, labels, False)

weights = ['3', '4']
sizes = ['20', '30']
fileNames = ['ldosBBH3D_L' + size + '_w' + w + '_E0_nMu1024_m1.1_diag' for w in weights for size in sizes]
labels = ['W = ' + w + ', L = ' + size for w in weights for size in sizes]
plots("ldosBBH3D_2D_E0_nMu1024_diag", fileNames, labels, False)
