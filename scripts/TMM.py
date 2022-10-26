import helper as hp
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

plt.style.use('science')

def plotTMM(name, label, ax):
    data = hp.readfile(name + ".dat")
    data = data[:, data[0,:].argsort()]

    ax.plot(data[0], data[1], label = label, linestyle='-')
def plotConstW(name, fileNames, labels, xlabel, show = True):
    fig, ax = plt.subplots()

    for i in range(0, len(fileNames)):
        plotTMM(fileNames[i], labels[i], ax)

    ax.set(xlabel = xlabel, ylabel = r'$\Lambda$')
    plt.xscale('log')
    plt.yscale('log')
    ax.legend(fontsize = 6, ncol = 2)
    ax.xaxis.set_minor_formatter(mticker.ScalarFormatter())
    #ax.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
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
    ax.margins(x = 0)
    ax.set_xlim([0,9])
    #hp.sotaiPhases(ax, 0.75)
    #ax.legend(fontsize = 6, ncol = 1, bbox_to_anchor=(0.7,0.6))
    ax.legend(fontsize = 6, ncol = 2)

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 300)
    fig.savefig(hp.plot_dir() + name + ".eps")
    if(show):
        plt.show()
    plt.close()


constWVals = ["-3", "-2.5", "-2", "-1.5", "-1", "-0.5", "0"]
constWNames = ["tmmSOTAI_E" + val + "_w4.6_d1_m1.1" for val in constWVals]
constWLabels= ["E = " + val for val in constWVals]

#plotConstW("tmmSOTAI_w4.6_d1_m1.1", constWNames, constWLabels, r'$L_x$', False)

constWVals = ["-3", "-2.5", "-2", "-1.5", "-1", "-0.5", "0"]
constWNames = ["tmmSOTAI_E" + val + "_w3.2_d1_m1.1" for val in constWVals]
constWLabels= ["E = " + val for val in constWVals]

#plotConstW("tmmSOTAI_w3.2_d1_m1.1", constWNames, constWLabels, r'$L_x$', False)

constLVals = ["20", "40", "60", "80", "100"]
constLNames = ["tmmSOTAI_E0_L" + val + "_d1_m1.1" for val in constLVals]
constLLabels = [r'$L_x = $ ' + val for val in constLVals]

#plotConstL("tmmSOTAI_E0_d1_m1.1_v1", constLNames, constLLabels, False)

constWVals = ["-3", "-2.5", "-2", "-1.5", "-1", "-0.5", "0"]
constWNames = ["tmmSOTAI_E" + val + "_w4.6_d0_m1.1" for val in constWVals]
constWLabels= ["E = " + val for val in constWVals]

#plotConstW("tmmSOTAI_w4.6_d0_m1.1", constWNames, constWLabels, r'$L_y$', False)

constWVals = ["-3", "-2.5", "-2", "-1.5", "-1", "-0.5", "0"]
constWNames = ["tmmSOTAI_E" + val + "_w3.2_d0_m1.1" for val in constWVals]
constWLabels= ["E = " + val for val in constWVals]

#plotConstW("tmmSOTAI_w3.2_d0_m1.1", constWNames, constWLabels, r'$L_y$', False)

constLVals = ["20"]
constLNames = ["tmmSOTAI_E0_L" + val + "_d0_m1.1" for val in constLVals]
constLLabels = [r'$L_y = $ ' + val for val in constLVals]

#plotConstL("tmmSOTAI_E0_d0_m1.1", constLNames, constLLabels, False)

constLVals = ["20", "40", "60", "80", "100"]
constLNames = ["tmmSOTAI_E0_L" + val + "_d1_m1.1" for val in constLVals]
constLLabels = [r'$L_x = $ ' + val for val in constLVals]

#plotConstL("tmmSOTAI_E0_d1_m1.1", constLNames, constLLabels, False)

constLVals = ["2", "4", "5", "6", "8", "10"]
constLNames = ["tmmBBH3D_E0_L" + val + "_d2_m1.1" for val in constLVals]
constLLabels = [r'$L_{x/y} = $ ' + val for val in constLVals]

plotConstL("tmmBBH3D_E0_d2_m1.1", constLNames, constLLabels, False)
