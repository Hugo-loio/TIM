import helper as hp
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

plt.style.use('science')

def probDens(ax, fname, label, diag: bool):
    data = hp.readfile(fname)
    x = data[0::4,0]
    y = data[1::4,0]
    z = data[2::4,0]
    d = np.sqrt(x*x + y*y + z*z)

    prob = np.average(data[3::4], axis = 1)
    prob_err = np.std(data[3::4], axis = 1)/np.sqrt(len(data[0]))

    if(diag):
        i = 0
        j = int(len(d)/2)
    else:
        i = int(len(d)/2)
        j = int(len(d))
    d = d[i:j]
    prob = prob[i:j]
    prob_err = prob_err[i:j]

    #ax.plot(d, dos, label = label)
    ax.errorbar(d, prob, yerr = prob_err, capsize = 2, linewidth = 0.5, capthick = 0.5, label = label)


def plots(name, fileNames, labels, show, diag):
    fig, ax = plt.subplots()
    ax.set(xlabel = r'$| \overrightarrow{n} |$', ylabel = r'$| \psi(\overrightarrow{n})|^2$')

    for i in range(len(fileNames)):
        probDens(ax, fileNames[i] + ".dat", labels[i], diag)

    #plt.yscale('log')
    #plt.xscale('log')
    #ax.margins(x = 0)
    ax.legend(fontsize = 6, ncol = 1)

    fig.savefig(hp.plot_dir() + name, dpi = 300, bbox_inches='tight')
    if(show):
        plt.show()

weights = ['3', '4']
sizes = ['20']
fileNames = ['probDensBBH3D_L' + size + '_w' + w + '_E0_m1.1' for w in weights for size in sizes]
labels = ['W = ' + w + ', L = ' + size for w in weights for size in sizes]
plots("probDensBBH3D_E0_edge", fileNames, labels, False, False)
plots("probDensBBH3D_E0_diag", fileNames, labels, False, True)
