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


def plotEnGap(name, fileName, show):
    fig, ax = plt.subplots()
    ax.set(xlabel = r'$W$', ylabel = r'Energy')

    data = hp.readfile(fileName + ".dat")
    w = data[0]
    gap = []
    gap_err = []
    for i in range(len(data[0])):
        gaps = [data[j+1,i] - data[j,i] for j in range(1, len(data) -1) if (data[j+1,i] > 0 and data[j,i] < 0) ]
        print(len(gaps))
        print(len(data[:,i]))
        gap.append(np.average(gaps))
        gap_err.append(np.std(gaps)/np.sqrt(len(gaps)))

    #plt.yscale('log')
    #plt.xscale('log')
    #ax.margins(x = 0)
    #ax.legend(fontsize = 6, ncol = 1)

    ax.errorbar(w, gap, yerr = gap_err, capsize = 2, linewidth = 0.5, capthick = 0.5)

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 300, bbox_inches='tight')
    if(show):
        plt.show()

sizes = ['20']
fileNames = ['spectrumE0BBH3D_L' + size + '_m1.1' for size in sizes]
for fileName in fileNames:
    plotEnGap(fileName + "Gap", fileName, False)
