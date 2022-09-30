import helper as hp
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

plt.style.use('science')

def plotTMM(name, label, ax):
    data = hp.readfile(name + ".dat")
    data = data[:, data[0,:].argsort()]

    ax.plot(data[0], data[1], label = label, linestyle='-')
def plotConstW(name, fileNames, labels, show = True):
    fig, ax = plt.subplots()

    for i in range(0, len(fileNames)):
        plotTMM(fileNames[i], labels[i], ax)

    ax.set(xlabel = r'$L_x$', ylabel = r'$\Lambda$')
    plt.xscale('log')
    plt.yscale('log')
    ax.legend()
    ax.xaxis.set_minor_formatter(mticker.ScalarFormatter())
    ax.yaxis.set_minor_formatter(mticker.ScalarFormatter())

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 300)
    fig.savefig(hp.plot_dir() + name + ".eps")
    if(show):
        plt.show()
    plt.close()

def plotConstL(name, fileNames, labels, show = True):
    fig, ax = plt.subplots()

    for i in range(0, len(fileNames)):
        plotTMM(fileNames[i], labels[i], ax)

    ax.set(xlabel = r'$W$', ylabel = r'$\Lambda$')
    plt.yscale('log')
    ax.legend()

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 300)
    fig.savefig(hp.plot_dir() + name + ".eps")
    if(show):
        plt.show()
    plt.close()


constW1Names = ["3.2", "0"]
constW1Names = ["tmmSOTAI_E" + name + "_w3.2_m1.1" for name in constW1Names]
constW1Labels= ["E = 3.2", "E = 0"]

#plotConstW("tmmSOTAI_w3.2_m1.1", constW1Names, constW1Labels, False)

constLVals = ["10", "16", "20", "30", "40", "60", "80", "100", "120"]
constLNames = ["tmmSOTAI_E0_L" + val + "_m1.1" for val in constLVals]
constLLabels = [r'$L_x = $ ' + val for val in constLVals]

plotConstL("tmmSOTAI_E0_m1.1", constLNames, constLLabels, True)
