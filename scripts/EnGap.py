import helper as hp
import matplotlib.pyplot as plt
import numpy as np

plt.style.use('science')

def plotEnGap(name, label, ax, detail):
    data = hp.readfile(name + ".dat")
    nPoints = int(len(data)/2)
    rho = np.average(data[2::2], axis = 0)
    #print(rho)
    #print(np.average(rho))

    if(detail):
        nPoints = len(data[2::2,0])
        errRho = np.std(data[2::2], axis = 0)/np.sqrt(nPoints)
        ax.errorbar(data[0], rho, label = label, yerr = errRho , capsize = 5, linestyle='-', linewidth=0.5)
    else:
        ax.plot(data[0], rho, label = label, linestyle='-', linewidth=0.5, alpha = 1)

def plot(name, fileNames, labels, detail = True, show = True):
    fig, ax = plt.subplots()

    for i in range(0, len(fileNames)):
        plotEnGap(fileNames[i], labels[i], ax, detail)

    ax.set(xlabel = r'$W$', ylabel = r'Gap')
    ax.legend(loc = 'upper right', fontsize = 7)

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 300)
    fig.savefig(hp.plot_dir() + name + ".eps")
    if(show):
        plt.show()
    plt.close()

size = ["100"]
mu = ["8192"]
rand = ["1"]
names = ["dosSOTAI_L" + size[i] + "_E0_nMu" + mu[i] + "_nR" + rand[i] + "_m1.1" for i in range(len(size))]
labels = ["L = " + size[i] + ", N = " + mu[i] + ", R = " + rand[i] for i in range(len(size))]

plot("EnGapSOTAI_intra1.1_test", names, labels, False, True)
