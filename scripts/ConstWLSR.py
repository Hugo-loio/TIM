import helper as hp
import matplotlib.pyplot as plt
import numpy as np

plt.style.use('science')

def plotLSR(name, label, ax, detail):
    data = hp.readfile(name + ".dat")
    data = data[:, data[0,:].argsort()]

    if(detail):
        ax.errorbar(data[0], np.average(data[1:],axis=0), label = label, yerr = np.std(data[1:], axis = 0)/np.sqrt(len(data[1:])) , capsize = 5, linestyle='-')
    else:
        ax.plot(data[0], np.average(data[1:],axis=0), label = label, linestyle='-')

def plot(name, fileNames, labels, detail = True, show = True):
    fig, ax = plt.subplots()

    for i in range(0, len(fileNames)):
        plotLSR(fileNames[i], labels[i], ax, detail)

    ax.set(xlabel = r'$L$', ylabel = r'LSR')
    ax.legend()

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 300)
    fig.savefig(hp.plot_dir() + name + ".eps")
    if(show):
        plt.show()
    plt.close()


weights = ["0.1"]
nstates = ["50"]
names = ["lsrSOTAI_w" + weights[i] + "_n" + nstates[i] + "_m1.1" for i in range(len(weights))]
labels = ["w = " + weights[i] + ", n = " + nstates[i] for i in range(len(weights))]
plot("ConstWLSR_SOTAI_intra1.1", names, labels, True, True)
