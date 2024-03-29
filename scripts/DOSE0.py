import helper as hp
import matplotlib.pyplot as plt
import numpy as np

plt.style.use('science')

def plotDOS(name, label, ax, detail):
    data = hp.readfile(name + ".dat")
    rho = np.average(data[1:], axis = 0)
    #print(rho)
    #print(np.average(rho))

    if(detail):
        nPoints = len(data[1:,0])
        errRho = np.std(data[1:], axis = 0)/np.sqrt(nPoints)
        ax.errorbar(data[0], rho, label = label, yerr = errRho , capsize = 5, linestyle='-', linewidth=0.5)
    else:
        ax.plot(data[0], rho, label = label, linestyle='-', linewidth=0.5, alpha = 1)

def plot(name, fileNames, labels, detail = True, show = True):
    fig, ax = plt.subplots()

    for i in range(0, len(fileNames)):
        plotDOS(fileNames[i], labels[i], ax, detail)

    ax.set(xlabel = r'$W$', ylabel = r'$\rho(0)$')
    ax.margins(x = 0)
    if(len(fileNames) > 1):
        #ax.legend(fontsize = 8, title = 'N')
        ax.legend(title = 'N', bbox_to_anchor = (0.6, 0.7))
    hp.sotaiPhases(ax, 0.85)
    #hp.totaiPhases(ax,0.8)

    #fig.set_size_inches(2.3,1.7)
    fig.set_size_inches(4,3)
    fig.savefig(hp.plot_dir() + name + ".png", bbox_inches = 'tight', dpi = 300)
    fig.savefig(hp.plot_dir() + name + ".pdf", bbox_inches = 'tight', pad_inches = 0.01)
    if(show):
        plt.show()
    plt.close()

size = ["100"]
mu = ["8192"]
rand = ["1"]
names = ["dosSOTAI_L" + size[i] + "_E0_nMu" + mu[i] + "_nR" + rand[i] + "_m1.1" for i in range(len(size))]
labels = ["L = " + size[i] + ", N = " + mu[i] + ", R = " + rand[i] for i in range(len(size))]

plot("DOSE0SOTAI_intra1.1", names, labels, False, False)

size = ["80", "80"]
mu = ["2048", "4096"]
rand = ["1", "1"]
names = ["dosBBH3D_L" + size[i] + "_E0_nMu" + mu[i] + "_nR" + rand[i] + "_m1.1" for i in range(len(size))]
labels = [mu[i]  for i in range(len(size))]

#plot("DOSE0BBH3D_intra1.1", names, labels, False, False)
