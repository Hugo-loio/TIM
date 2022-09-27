import helper as hp
import matplotlib.pyplot as plt
import numpy as np

plt.style.use('science')

def plotDOS(name, label, ax, detail):
    data = hp.readfile(name + ".dat")
    nPoints = int(len(data)/2)
    en = data[0::2,0]
    rho = np.average(data[1::2,:], axis = 1)
    #print(rho)
    #print(np.average(rho))

    if(detail):
        errRho = np.std(data[nPoints:,:], axis = 1)/np.sqrt(len(data[0]))
        ax.errorbar(en, rho, label = label, yerr = errRho , capsize = 5, linestyle='-')
    else:
        ax.plot(en, rho, label = label, linestyle='-', linewidth=0.5, alpha = 1)

def plot(name, fileNames, labels, detail = True, show = True):
    fig, ax = plt.subplots()

    for i in range(0, len(fileNames)):
        plotDOS(fileNames[i], labels[i], ax, detail)

    ax.set(xlabel = r'$E$', ylabel = r'$\rho(E)$')
    ax.legend(loc = 'upper right', fontsize = 7)

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 300)
    fig.savefig(hp.plot_dir() + name + ".eps")
    if(show):
        plt.show()
    plt.close()

def plotZoom(name, fileNames, labels, detail = True, show = True):
    fig, ax = plt.subplots()

    for i in range(0, len(fileNames)):
        plotDOS(fileNames[i], labels[i], ax, detail)

    ax.set(xlabel = r'$E$', ylabel = r'$\rho(E)$')
    ax.legend(loc = 'upper right', fontsize = 7)
    plt.xlim([-0.5,0.5])

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 300)
    fig.savefig(hp.plot_dir() + name + ".eps")
    if(show):
        plt.show()
    plt.close()


size = ["100", "100", "200", "200", "200", "200", "200"]
weight = ["3.2", "3.2", "3.2", "3.2", "3.2", "3.2", "3.2"]
mu = ["1000", "1000", "1000", "1000", "100", "100", "10000"]
rand = ["1", "5", "5", "1", "1", "5", "1", "1"]
names = ["dosSOTAI_L" + size[i] + "_w" + weight[i] + "_nMu" + mu[i] + "_nR" + rand[i] + "_m1.1" for i in range(len(size))]
labels = ["L = " + size[i] + ", W = " + weight[i] + ", N = " + mu[i] + ", R = " + rand[i] for i in range(len(size))]

#plot("ConstWDosSOTAI_intra1.1_test", names, labels, False, False)

weight = ["2.4", "2.8", "3.2", "3.6", "4", "9"]
names = ["dosSOTAI_L200" + "_w" + weight[i] + "_nMu1000_nR1_m1.1" for i in range(len(weight))]
labels = ["W = " + weight[i] for i in range(len(weight))]
#plot("ConstWDosSOTAI_intra1.1_v1", names, labels, False, False)
#plotZoom("ConstWDosSOTAI_intra1.1_v1_zoom", names, labels, False, False)

names = ["dosSOTAI_L100" + "_w" + weight[i] + "_nMu8192_nR1_m1.1" for i in range(len(weight))]
plot("ConstWDosSOTAI_intra1.1_v2", names, labels, False, False)
plotZoom("ConstWDosSOTAI_intra1.1_v2_zoom", names, labels, False, False)

size = ["50", "50"]
weight = ["1", "2"]
mu = ["100", "100"]
rand = ["1", "1"]
names = ["dosSOTAI_L" + size[i] + "_w" + weight[i] + "_nMu" + mu[i] + "_nR" + rand[i] + "_m1.1" for i in range(len(size))]
labels = ["L = " + size[i] + ", W = " + weight[i] + ", N = " + mu[i] + ", R = " + rand[i] for i in range(len(size))]

#plot("ConstWDosSOTAI_intra1.1_test", names, labels, False, False)
#plotZoom("ConstWDosSOTAI_intra1.1_v1_zoom", names, labels, False, False)


