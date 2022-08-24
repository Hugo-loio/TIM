import helper as hp
import matplotlib.pyplot as plt
import numpy as np

def plotTMM(name, label, ax):
    data = hp.readfile(name + ".dat")
    data = data[:, data[0,:].argsort()]

    ax.plot(data[0], data[1], label = label, linestyle='-')

def plotConstW(name, fileNames, labels, show = True):
    fig, ax = plt.subplots()

    for i in range(0, len(fileNames)):
        plotTMM(fileNames[i], labels[i], ax)

    ax.set(xlabel = r'$L_x$', ylabel = r'$\Lambda$')
    ax.legend()

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 200)
    if(show):
        plt.show()
    plt.close()

def plotConstL(name, fileNames, labels, show = True):
    fig, ax = plt.subplots()

    for i in range(0, len(fileNames)):
        plotTMM(fileNames[i], labels[i], ax)

    ax.set(xlabel = r'$W$', ylabel = r'$\Lambda$')
    ax.legend()

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 200)
    if(show):
        plt.show()
    plt.close()


constW1Names = ["0"]
constW1Names = ["tmmSOTAI_E" + name + "_w3.2_m1.1" for name in constW1Names]
constW1Labels= ["E = 0"]

plotConstW("tmmSOTAI_w3.2_m1.1", constW1Names, constW1Labels, False)

constLVals = ["20"]
constLNames = ["tmmSOTAI_E0_L" + val + "_m1.1" for val in constLVals]
constLLabels = ["L = 20"]

plotConstL("tmmSOTAI_E0_m1.1", constLNames, constLLabels, True)
