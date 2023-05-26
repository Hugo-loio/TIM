import helper as hp
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

plt.style.use('science')

def plotTMM(name, label, ax, mode = 0):
    data = hp.readfile(name + ".dat")
    data = data[:, data[0,:].argsort()]

    if(mode == 0):
        ax.plot(data[0], data[1], label = label, linestyle='-')
    elif(mode == 1):
        ax.plot(data[0][1:], data[1][1:], label = label, linestyle='-', marker='.')
        #ax.plot(data[0][1:], data[1][1:], label = label, linestyle='-', marker='.', lw = 1, markersize = 2)

def plotConstW(name, fileNames, labels, xlabel, show = True, mode = 0, model = 0):
    fig, ax = plt.subplots()

    for i in range(0, len(fileNames)):
        plotTMM(fileNames[i], labels[i], ax, mode)

    ax.set(xlabel = xlabel, ylabel = r'$\Lambda$')
    plt.xscale('log')
    plt.yscale('log')
    #ax.legend(fontsize = 3.4, ncol = 2, title = r'$E$', title_fontsize = 7)
    if(model == 2):
        ax.legend(fontsize = 8, ncol = 2, title = r'$E$', title_fontsize = 8)
        ymin, ymax = ax.get_ylim()
        #ax.set_ylim([ymin, ymax*1.2])
        ax.yaxis.labelpad = 0
        ax.tick_params(axis='y', which='major', pad=0)
        #ax.text(6.4,4, "Phase III")
        ax.text(7.5,4, "Phase III")
        #ax.text(4,11, "Phase II")
    elif(model == 1):
        ax.legend(fontsize = 9, ncol = 2, title = r'$E$')
        #ax.yaxis.labelpad = 0
        #ax.tick_params(axis='y', which='major', pad=0)
        #ax.text(50,1, "Phase III")
        ax.text(50,1, "Phase IV")

    ax.yaxis.set_minor_formatter(ticker.ScalarFormatter()) 
    ax.yaxis.set_minor_locator(MultipleLocator(0.4))
    #ax.yaxis.set_minor_formatter(ticker.NullFormatter()) 
    #ax.yaxis.set_major_formatter(ticker.ScalarFormatter()) 

    #ax.xaxis.set_major_formatter(ticker.ScalarFormatter()) 
    ax.xaxis.set_minor_formatter(ticker.ScalarFormatter()) 
    ax.xaxis.set_minor_locator(MultipleLocator(20))

    #fig.set_size_inches(1.7,1.3)
    fig.set_size_inches(4,3)
    fig.savefig(hp.plot_dir() + name + ".png", bbox_inches = 'tight', dpi = 300)
    fig.savefig(hp.plot_dir() + name + ".pdf", bbox_inches = 'tight', pad_inches = 0.01)
    if(show):
        plt.show()
    plt.close()

def plotConstL(name, fileNames, labels, show = True, model = 0):
    fig, ax = plt.subplots()

    for i in range(0, len(fileNames)):
        plotTMM(fileNames[i], labels[i], ax)

    ax.set(xlabel = r'$W$', ylabel = r'$\Lambda$')
    plt.yscale('log')
    ax.margins(x = 0)
    if(model == 1):
        hp.sotaiPhases(ax, 0.75)
        ax.legend(fontsize = 8, ncol = 1, bbox_to_anchor=(0.68,0.5))
    elif(model == 2):
        ymin, ymax = ax.get_ylim()
        ax.set_xlim([0,9])
        ax.set_ylim([ymin, ymax*1.3])
        #hp.totaiPhases(ax, 0.5)
        hp.totaiPhases(ax, 0.7)
        #ax.legend(fontsize = 4, ncol = 2, title = r'$L$', title_fontsize = 7)
        ax.legend(fontsize = 8, ncol = 2, title = r'$L$', title_fontsize = 9, bbox_to_anchor = (0.7, 0.15), loc = 'center')
        ax.tick_params(axis='y', which='major', pad=0)
    elif(model == 3):
        ax.set_xlim([18,34])
        ax.set_ylim([0.2,3])
        #ax.legend(fontsize = 5, ncol = 2, title = r'$L$', title_fontsize = 8, bbox_to_anchor = (0.7,0.65), loc = 'center')
        ax.legend(fontsize = 8, ncol = 2, title = r'$L$', title_fontsize = 9, bbox_to_anchor = (0.7, 0.65), loc = 'center')
        ax.tick_params(axis='y', which='major', pad=0)
        ax.yaxis.set_minor_formatter(ticker.NullFormatter()) 
        #ax.yaxis.set_minor_formatter(ticker.ScalarFormatter()) 
        ax.axvline(24, color='black', linestyle='dashed', linewidth=0.8)
        ax.text((18+24)/2, 2.2, r'III', ha = 'center', fontsize=8)
        ax.text((24+34)/2, 2.2, r'IV', ha = 'center', fontsize=8)
        ax.tick_params(axis='y', which='major', pad=0)
        #ax.legend(fontsize = 6, ncol = 2)
    else:
        ax.legend(fontsize = 6, ncol = 2)

    ax.yaxis.labelpad = 0

    #fig.set_size_inches(1.7,1.3)
    fig.set_size_inches(4,3)
    #fig.set_size_inches(2.3,3)
    fig.savefig(hp.plot_dir() + name + ".png", bbox_inches = 'tight', dpi = 300)
    fig.savefig(hp.plot_dir() + name + ".pdf", bbox_inches = 'tight', pad_inches = 0.01)
    if(show):
        plt.show()
    plt.close()


constWVals = ["-3", "-2.5", "-2", "-1.5", "-1", "-0.5", "0"]
constWNames = ["tmmSOTAI_E" + val + "_w4.6_d1_m1.1" for val in constWVals]
constWLabels= [val for val in constWVals]

plotConstW("tmmSOTAI_w4.6_d1_m1.1", constWNames, constWLabels, r'$L_x$', False, 0 , 1)

constWVals = ["-3", "-2.5", "-2", "-1.5", "-1", "-0.5", "0"]
constWNames = ["tmmSOTAI_E" + val + "_w3.2_d1_m1.1" for val in constWVals]
constWLabels= [val for val in constWVals]

plotConstW("tmmSOTAI_w3.2_d1_m1.1", constWNames, constWLabels, r'$L_x$', False, 0, 1)

constLVals = ["20", "40", "60", "80", "100"]
constLNames = ["tmmSOTAI_E0_L" + val + "_d1_m1.1" for val in constLVals]
constLLabels = [r'$L_x = $ ' + val for val in constLVals]

#plotConstL("tmmSOTAI_E0_d1_m1.1_v1", constLNames, constLLabels, False, 1)

constWVals = ["-3", "-2.5", "-2", "-1.5", "-1", "-0.5", "0"]
constWNames = ["tmmSOTAI_E" + val + "_w4.6_d0_m1.1" for val in constWVals]
constWLabels= ["E = " + val for val in constWVals]

#plotConstW("tmmSOTAI_w4.6_d0_m1.1", constWNames, constWLabels, r'$L_y$', False)

constWVals = ["-3", "-2.5", "-2", "-1.5", "-1", "-0.5", "0"]
constWNames = ["tmmSOTAI_E" + val + "_w3.2_d0_m1.1" for val in constWVals]
constWLabels= ["E = " + val for val in constWVals]

#plotConstW("tmmSOTAI_w3.2_d0_m1.1", constWNames, constWLabels, r'$L_y$', False)

constWVals = ["-3", "-2.5", "-2", "-1.5", "-1", "-0.5", "0"]
constWNames = ["tmmBBH3D_E" + val + "_w3.4_d2_m1.1" for val in constWVals]
constWLabels= [val for val in constWVals]

#plotConstW("tmmBBH3D_w3.4_d2_m1.1", constWNames, constWLabels, r'$L$', False, 1)

constWVals = ["-3", "-2.5", "-2", "-1.5", "-1", "-0.5", "0"]
constWNames = ["tmmBBH3D_E" + val + "_w5_d2_m1.1" for val in constWVals]
constWLabels= [val for val in constWVals]

#plotConstW("tmmBBH3D_w5_d2_m1.1", constWNames, constWLabels, r'$L$', False, 1)

constWVals = ["-3", "-2.5", "-2", "-1.5", "-1", "-0.5", "0"]
constWNames = ["tmmBBH3D_E" + val + "_w50_d2_m1.1" for val in constWVals]
constWLabels= ["E = " + val for val in constWVals]

#plotConstW("tmmBBH3D_w50_d2_m1.1", constWNames, constWLabels, r'$L_{x/y}$', False, 1)

constLVals = ["20"]
constLNames = ["tmmSOTAI_E0_L" + val + "_d0_m1.1" for val in constLVals]
constLLabels = [r'$L_y = $ ' + val for val in constLVals]

#plotConstL("tmmSOTAI_E0_d0_m1.1", constLNames, constLLabels, False)

constLVals = ["20", "40", "60", "80", "100"]
constLNames = ["tmmSOTAI_E0_L" + val + "_d1_m1.1" for val in constLVals]
constLLabels = [r'$L_x = $ ' + val for val in constLVals]

#plotConstL("tmmSOTAI_E0_d1_m1.1", constLNames, constLLabels, False, 1)

constLVals = ["2", "4", "5", "6", "8", "10"]
constLNames = ["tmmBBH3D_E0_L" + val + "_d2_m1.1" for val in constLVals]
constLLabels = [val for val in constLVals]

#plotConstL("tmmBBH3D_E0_d2_m1.1", constLNames, constLLabels, False, 2)

dVals = ["0", "1", "2"]
constLNames = ["tmmBBH3D_E0_L4_d" + val + "_m1.1" for val in dVals]
constLLabels = ["x", "y", "z"]

#plotConstL("tmmBBH3D_E0_dComp_m1.1", constLNames, constLLabels, False)

constLVals = ["4", "6", "8", "10", "12"]
constLNames = ["tmmBBH3D_E0_L" + val + "_d2_m1.1_cross2" for val in constLVals]
constLLabels = [val for val in constLVals]

#plotConstL("tmmBBH3D_E0_d2_m1.1_cross2", constLNames, constLLabels, False, 3)
