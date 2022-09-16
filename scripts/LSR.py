import helper as hp
import matplotlib.pyplot as plt
import numpy as np

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

    ax.set(xlabel = r'$W$', ylabel = r'LSR')
    ax.legend()

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 300)
    fig.savefig(hp.plot_dir() + name + ".eps")
    if(show):
        plt.show()
    plt.close()


names = ["L50_n10", "L50_n10", "L100_n10"]
names = ["lsrSOTAI_" + name + "_m1.1" for name in names]
names[1] += "_open"
labels= ["L = 50, n = 10", "L = 50, n = 10, OBC", "L = 100, n = 10"]

#plot("lsrSOTAI_intra1.1_test", names, labels, True, True)

sizes = ["50", "100"]
nstates = ["50", "50"]
names = ["lsrSOTAI_L" + sizes[i] + "_n" + nstates[i] + "_m1.1" for i in range(len(sizes))]
labels = ["L = " + sizes[i] + ", n = " + nstates[i] for i in range(len(sizes))]
plot("LSR_SOTAI_intra1.1", names, labels, True, True)
