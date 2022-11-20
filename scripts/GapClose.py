import helper as hp
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

plt.style.use('science')

def plot(title, fileName, specs, errors = False, show = False):
    data = hp.readfile(fileName + ".dat")
    step = specs['nStSamp']*2 + 1

    fig, ax = plt.subplots()
    for w in specs['wVals']:
        y = []
        yerr = []
        for en in specs['enGapSizes']:
            dataEn = data[:,np.where(data[0] == i)[0]]
            dataW = dataEn[:np.where(data[1] == w)[0]]
            y.append(

    if(len(specs['enGapSizes']) > 0):
        fig, ax = plt.subplots()
        for i in specs['enGapSizes']:
            dataPlot = data[:,np.where(data[0] == i)[0]] 
            start = int(2 + 2*specs['nStSamp'])
            label = r'$L = $' + str(i)
            plotCurve(ax, dataPlot[1], dataPlot[start::step], label, errors)

        ax.set(xlabel = r'$W$', ylabel = r'Gap')
        ax.margins(x = 0)
        if(len(specs['enGapSizes']) > 1): 
            ax.legend(fontsize = 6)
        if(specs['phases'] == 0):
            hp.sotaiPhases(ax)

        fig.savefig(hp.plot_dir() + title + "Gap.png", dpi = 300)
        fig.savefig(hp.plot_dir() + title + "Gap.eps")
        if(show):
            plt.show()
        plt.close()

specs = {
        'nStSamp' : 5,
        'enGapSizes' : [6,8,10,12,14,16,18,20],
        'wVals' : [3.69, 3.78]
        }

plot("GapCloseBBH3D", "locBBH3D_m1.1", specs, False, True)
#plot("LocBBH3D_intra1.1", "locBBH3D_intra1.1", specs)
