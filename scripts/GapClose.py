import helper as hp
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

plt.style.use('science')

def plot(title, fileName, specs, errors = False, show = False):
    data = hp.readfile(fileName + ".dat")
    step = specs['nStSamp']*2 + 1
    start = int(2 + 2*specs['nStSamp'])

    fig, ax = plt.subplots()
    for w in specs['wVals']:
        y = []
        yerr = []
        for en in specs['enGapSizes']:
            dataEn = data[:,np.where(data[0] == en)[0]]
            dataEnW = dataEn[:,np.where(dataEn[1] == w)[0]][start::step]
            y.append(np.average(dataEnW, axis = 0)[0])
            yerr.append(np.std(dataEnW, axis = 0)[0]/np.sqrt(len(dataEnW)))

        ax.errorbar(specs['enGapSizes'], y, yerr = yerr, label = "W = " + str(w), capsize = 5, linestyle='-', linewidth = 0.5)

    ax.set(xlabel = r'$L$', ylabel = r'Gap')
    ax.margins(x = 0)
    ax.legend(fontsize = 6)
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.show()

specs = {
        'nStSamp' : 5,
        'enGapSizes' : [6,8,10,12,14,16,18,20],
        'wVals' : [3.69, 3.78]
        }

plot("GapCloseBBH3D", "locBBH3D_m1.1", specs, False, True)
#plot("LocBBH3D_intra1.1", "locBBH3D_intra1.1", specs)
