import helper as hp
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

plt.style.use('science')

def plotEnGap(name, fileName, show):
    data = hp.readfile(fileName + ".dat")
    w = data[0]
    gap = []
    gapErr = []
    for i in range(len(data[0])):
        gaps = [data[j+1,i] - data[j,i] for j in range(1, len(data) -1) if (data[j+1,i] > 0 and data[j,i] < 0) ]
        gap.append(np.average(gaps))
        gapErr.append(np.std(gaps)/np.sqrt(len(gaps)))
    mls = []
    mlsErr = []
    nums = [10,20,30,40,50]
    nmax = 50
    nsamp = len(data)/nmax
    midl = nmax/2-1
    for i in range(len(data[0])):
        ls = [data[j+1,i] - data[j,i] for j in range(1, len(data) -1) if (data[j+1,i]*data[j,i] > 0) ]
        mls.append([])
        mlsErr.append([])
        for n in nums:
            n2 =n/2-1
            lsN = np.array([ls[int(midl*i-n2):int(midl*i+n2)] for i in range(1, int(nsamp+1))]).flatten()
            mls[i].append(np.average(lsN))
            mlsErr[i].append(np.std(lsN)/np.sqrt(len(lsN)))
    mls = np.array(mls)
    mlsErr = np.array(mlsErr)

    fig, ax = plt.subplots()
    ax.set(xlabel = r'$W$', ylabel = r'Energy')

    plt.yscale('log')
    #plt.xscale('log')
    #ax.margins(x = 0)

    #ax.errorbar(w, gap, yerr = gapErr, capsize = 2, linewidth = 0.5, capthick = 0.5, label = "Gap")
    ax.errorbar(w, gap, capsize = 2, linewidth = 0.5, capthick = 0.5, label = "Gap")
    #ax.errorbar(w, mls, yerr = mlsErr, capsize = 2, linewidth = 0.5, capthick = 0.5, label = "MLS")
    for i,n in enumerate(nums):
        ax.errorbar(w, mls[:,i], capsize = 2, linewidth = 0.5, capthick = 0.5, label = "MLS, n = " + str(n))
    #ax.errorbar(w, mls, capsize = 2, linewidth = 0.5, capthick = 0.5, label = "MLS with gap")
    ax.legend(fontsize = 10, ncol = 1)
    #ax.legend()

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 300, bbox_inches='tight')
    if(show):
        plt.show()

sizes = ['20']
fileNames = ['spectrumE0BBH3D_L' + size + '_m1.1' for size in sizes]
for fileName in fileNames:
    plotEnGap(fileName + "Gap", fileName, False)
