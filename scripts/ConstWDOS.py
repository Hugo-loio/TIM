import helper as hp
import matplotlib.pyplot as plt
import numpy as np

plt.style.use('science')

def plotDOS(name, label, ax, detail):
    data = hp.readfile(name + ".dat")
    nPoints = int(len(data)/2)
    en = data[0:nPoints,0]
    rho = np.average(data[nPoints:,:], axis = 1)
    #print(rho)
    print(np.average(rho))

    if(detail):
        errRho = np.std(data[nPoints:,:], axis = 1)/np.sqrt(len(data[0]))
        ax.errorbar(en, rho, label = label, yerr = errRho , capsize = 5, linestyle='-')
    else:
        ax.plot(en, rho, label = label, linestyle='-')

def plot(name, fileNames, labels, detail = True, show = True):
    fig, ax = plt.subplots()

    for i in range(0, len(fileNames)):
        plotDOS(fileNames[i], labels[i], ax, detail)

    ax.set(xlabel = r'$E$', ylabel = r'$\rho(E)$')
    ax.legend()

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 200)
    if(show):
        plt.show()
    plt.close()


size = ["100"]
weight = ["3.2"]
mu = ["100"]
rand = ["5"]
names = ["dosSOTAI_L" + size[i] + "_w" + weight[i] + "_nMu" + mu[i] + "_nR" + rand[i] + "_m1.1" for i in range(len(size))]
labels = ["L = " + size[i] + ", W = " + weight[i] + ", N = " + mu[i] + ", R = " + rand[i] for i in range(len(size))]

plot("ConstWDosSOTAI_intra1.1", names[0:1], labels[0:1], False, True)
