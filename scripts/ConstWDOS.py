import helper as hp
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

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

    colormap = plt.cm.turbo
    colors = [colormap(i) for i in np.linspace(0, 1, len(fileNames))]
    ax.set_prop_cycle('color', colors)

    for i in range(0, len(fileNames)):
        plotDOS(fileNames[i], labels[i], ax, detail)

    ax.set(xlabel = r'$E$', ylabel = r'$\rho(E)$')
    ax.legend(loc = 'upper right', fontsize = 5, title = r'$W$', title_fontsize = 7)

    fig.set_size_inches(2.3,1.7)
    fig.savefig(hp.plot_dir() + name + ".png", bbox_inches = 'tight', dpi = 300)
    fig.savefig(hp.plot_dir() + name + ".pdf", bbox_inches = 'tight', pad_inches = 0.01)
    if(show):
        plt.show()
    plt.close()

def plotZoom(name, fileNames, labels, detail = True, show = True, ylim = 0):
    fig, ax = plt.subplots()

    colormap = plt.cm.turbo
    colors = [colormap(i) for i in np.linspace(0, 1, len(fileNames))]
    ax.set_prop_cycle('color', colors)

    for i in range(0, len(fileNames)):
        plotDOS(fileNames[i], labels[i], ax, detail)

    ax.set(xlabel = r'$E$', ylabel = r'$\rho(E)$')
    ax.legend(loc = 'upper right', fontsize = 5, title = r'$W$', title_fontsize = 7, ncols = 2)

    plt.xlim([-0.5,0.5])
    if(ylim != 0):
        plt.ylim([-0.05*ylim, ylim])

    fig.set_size_inches(2.3,1.7)
    fig.savefig(hp.plot_dir() + name + ".png", bbox_inches = 'tight', dpi = 300)
    fig.savefig(hp.plot_dir() + name + ".pdf", bbox_inches = 'tight', pad_inches = 0.01)
    if(show):
        plt.show()
    plt.close()

def linFunc(x, a, b):
    return a*x + b

def plotFit(name, fileNames, labels, detail = True, show = True, ylim = 0):
    fig, ax = plt.subplots()

    colormap = plt.cm.turbo
    colors = [colormap(i) for i in np.linspace(0, 1, len(fileNames))]
    ax.set_prop_cycle('color', colors)

    for i in range(0, len(fileNames)):
        data = hp.readfile(fileNames[i] + ".dat")
        nPoints = int(len(data)/2)
        en = data[0::2,0]
        indices = np.argwhere(np.where((en > 0.04) & (en < 0.65), en, 0))
        en = en[indices]
        rho = data[1::2,0][indices]

        ax.plot(en, rho, label = labels[i], linestyle='-', linewidth=0.5, alpha = 1)
        
        xFit = np.array(np.log(en)[:,0])
        yFit = np.array(np.log(rho)[:,0])
        popt, pcov = optimize.curve_fit(linFunc, xFit, yFit)
        err = np.sqrt(np.diag(pcov))
        print(popt, err)

    ax.set(xlabel = r'$E$', ylabel = r'$\rho(E)$')
    ax.legend(fontsize = 10, title = r'$W$', title_fontsize = 12, ncols = 1)
    ax.set_xscale('log')
    ax.set_yscale('log')

    #plt.xlim([-0.5,0.5])
    #if(ylim != 0):
    #    plt.ylim([-0.05*ylim, ylim])

    #fig.set_size_inches(2.3,1.7)
    fig.savefig(hp.plot_dir() + name + ".png", bbox_inches = 'tight', dpi = 300)
    fig.savefig(hp.plot_dir() + name + ".pdf", bbox_inches = 'tight', pad_inches = 0.01)
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
#plot("ConstWDosSOTAI_intra1.1_v2", names, labels, False, False)
#plotZoom("ConstWDosSOTAI_intra1.1_v2_zoom", names, labels, False, False)

size = ["50", "50"]
weight = ["1", "2"]
mu = ["100", "100"]
rand = ["1", "1"]
names = ["dosSOTAI_L" + size[i] + "_w" + weight[i] + "_nMu" + mu[i] + "_nR" + rand[i] + "_m1.1" for i in range(len(size))]
labels = ["L = " + size[i] + ", W = " + weight[i] + ", N = " + mu[i] + ", R = " + rand[i] for i in range(len(size))]

#plot("ConstWDosSOTAI_intra1.1_test", names, labels, False, False)
#plotZoom("ConstWDosSOTAI_intra1.1_v1_zoom", names, labels, False, False)

weight = ["2", "2.8", "3", "3.2", "3.4", "3.6", "4", "9"]
names = ["dosBBH3D_L80_w" + weight[i] + "_nMu2048_nR1_intra1.1" for i in range(len(weight))]
labels = ["W = " + weight[i] for i in range(len(weight))]

#plot("ConstWDosBBH3D_intra1.1_L80_nMu2048_nR1", names, labels, False, False)
#plotZoom("ConstWDosBBH3D_intra1.1_L80_nMu2048_nR1_zoom", names, labels, False, False)

weight = ["2", "2.5", "2.6", "3.2", "3.4", "3.6", "4", "9"]
names = ["dosBBH3D_L80_w" + weight[i] + "_nMu4096_nR1_intra1.1" for i in range(len(weight))]
labels = [weight[i] for i in range(len(weight))]

plot("ConstWDosBBH3D_intra1.1_L80_nMu4096_nR1", names, labels, False, False)
plotZoom("ConstWDosBBH3D_intra1.1_L80_nMu4096_nR1_zoom", names, labels, False, False)

weight = ["2.5", "2.6"]
names = ["dosBBH3D_L80_w" + weight[i] + "_nMu4096_nR1_intra1.1" for i in range(len(weight))]
labels = [weight[i] for i in range(len(weight))]

#plot("ConstWDosBBH3D_intra1.1_L80_nMu4096_nR1", names, labels, False, False)
#plotFit("ConstWDosBBH3D_intra1.1_L80_nMu4096_nR1_fit", names, labels, False, False)

weight = ["2.5", "2.55", "2.6"]
names = ["dosBBH3D_L80_w" + weight[i] + "_nMu8192_nR1_intra1.1" for i in range(len(weight))]
labels = ["W = " + weight[i] for i in range(len(weight))]

#plot("ConstWDosBBH3D_intra1.1_L80_nMu8192_nR1", names, labels, False, False)
#plotZoom("ConstWDosBBH3D_intra1.1_L80_nMu8192_nR1_zoom", names, labels, False, False)

intra = ["0.5","0.9", "1", "1.1", "2"]
names = ["dosBBH3D_L80_intra" + i + "_w0_nMu4096_nR1" for i in intra]
labels = ["$\gamma$ = " + i for i in intra]

#plot("CleanDosBBH3D_L80_nMu4096_nR1", names, labels, False, False)
#plotZoom("CleanDosBBH3D_L80_nMu4096_nR1_zoom", names, labels, False, False, 0.03)

intra = ["0.5","0.9", "1", "1.1", "2"]
names = ["dosBBH3D_L80_intra" + i + "_w0_nMu2048_nR1" for i in intra]
labels = ["$\gamma$ = " + i for i in intra]

#plot("CleanDosBBH3D_L80_nMu2048_nR1", names, labels, False, False)
#plotZoom("CleanDosBBH3D_L80_nMu2048_nR1_zoom", names, labels, False, False, 0.03)

intra = ["0.5","0.9", "1", "1.1", "2"]
names = ["dosBBH3D_L80_intra" + i + "_w0_nMu1024_nR1" for i in intra]
labels = ["$\gamma$ = " + i for i in intra]

#plot("CleanDosBBH3D_L80_nMu1024_nR1", names, labels, False, False)
#plotZoom("CleanDosBBH3D_L80_nMu1024_nR1_zoom", names, labels, False, False, 0.03)
