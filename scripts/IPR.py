import helper as hp
import matplotlib.pyplot as plt
import numpy as np

def plotIPR(name, label, ax, detail):
    data = hp.readfile(name + ".dat")
    data = data[:, data[0,:].argsort()]

    if(detail):
        ax.errorbar(data[0], np.average(data,axis=0), label = label, yerr = np.std(data[1:], axis = 0)/np.sqrt(len(data[1:])) , capsize = 5, linestyle='-')
    else:
        ax.plot(data[0], np.average(data,axis=0), label = label, yerr = np.std(data[1:], axis = 0)/np.sqrt(len(data[1:])) , capsize = 5, linestyle='-')

def plot(name, fileNames, labels, detail = True, show = True):
    fig, ax = plt.subplots()

    for i in range(0, len(fileNames)):
        plotIPR(fileNames[i], labels[i], ax, detail)

    ax.set(xlabel = r'$W$', ylabel = r'IPR')
    ax.legend()

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 200)
    if(show):
        plt.show()
    plt.close()


names = ["L20_n10", "L50_n10", "L100_n10"]
names = ["iprSOTAI_" + name + "_m1.1" for name in names]
labels= ["L = 20, n = 10", "L = 50, n = 10", "L = 100, n = 10"]

plot("iprSOTAI_intra1.1", names, labels, True, True)
