import helper as hp
import matplotlib.pyplot as plt
import numpy as np

plt.style.use('science')

def plotEMax(name, label, ax):
    data = hp.readfile(name + ".dat")

    ax.plot(data[0], data[1], label = label, linestyle='-')

def plotSOTAI(name, fileNames, labels, show = True):
    fig, ax = plt.subplots()

    for i in range(0, len(fileNames)):
        plotEMax(fileNames[i], labels[i], ax)

    x = np.linspace(0, 9, 100)
    b = 4
    m = (11-4)/9
    f = m*x + b
    ax.plot(x, f, label = r'Safe limit')

    ax.set(xlabel = r'$W$', ylabel = r'E Max')
    ax.legend(fontsize = 8)

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 200)
    if(show):
        plt.show()
    plt.close()

def plotBBH3D(name, fileNames, labels, show = True):
    fig, ax = plt.subplots()

    for i in range(0, len(fileNames)):
        plotEMax(fileNames[i], labels[i], ax)

    x = np.linspace(0, 9, 100)
    b = 5
    m = 1
    f = m*x + b
    ax.plot(x, f, label = r'Safe limit')

    ax.set(xlabel = r'$W$', ylabel = r'$E_{\textnormal{max}}$')
    ax.legend(fontsize = 8)

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 300)
    if(show):
        plt.show()
    plt.close()


size = ["10", "10", "15", "20"]
nHam = ["100", "50", "50", "50"]
names = ["eMaxSOTAI_L" + size[i] + "_nHam" + nHam[i]  + "_m1.1" for i in range(len(size))]
labels = ["L = " + size[i] + " , n = " + nHam[i] for i in range(len(size))]

#plotSOTAI("EMaxSOTAI_intra1.1", names, labels, False)

size = ["5", "5", "8", "6", "8", "10"]
nHam = ["50", "100", "50", "50", "50", "50"]
names = ["eMaxBBH3D_L" + size[i] + "_nHam" + nHam[i]  + "_m1.1" for i in range(len(size))]
labels = ["L = " + size[i] + " , n = " + nHam[i] for i in range(len(size))]
plotBBH3D("EMaxBBH3D_intra1.1", names, labels, False)
