import helper as hp
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

plt.style.use('science')

def fitLin(x, a, b):
    return a*x + b

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

        x = 1/np.array(specs['enGapSizes'])
        sort = np.argsort(x)
        x = x[sort]
        y = np.array(y)[sort]
        yerr = np.array(yerr)[sort]

        y0 = []
        for i in range(3):
            x_temp = x[0:5-i]
            y_temp = y[0:5-i]
            #yerr_temp = yerr[0:5-i]
            #popt, pcov = optimize.curve_fit(fitLin, x_temp, y_temp, sigma = yerr_temp, absolute_sigma = True)
            popt, pcov = optimize.curve_fit(fitLin, x_temp, y_temp)
            y0.append(popt[1])
        
        y0_avg = np.average(y0)
        #print(y0_avg)
        y0_err = np.amax(abs(y0-y0_avg))
        x = np.insert(x,0,0)
        y = np.insert(y,0,y0_avg)
        yerr = np.insert(yerr,0,y0_err)

        ax.errorbar(x, y, yerr = yerr, label = "W = " + str(w), capsize = 2, linestyle='-', linewidth = 0.5, capthick = 0.5)

    ax.set(xlabel = r'$\frac{1}{L}$', ylabel = r'Gap')
    #ax.margins(x = 0)
    #ax.set_xlim([0, ax.get_xlim()[1]])
    ax.legend(fontsize = 6)
    #ax.set_xscale('log')
    #ax.set_yscale('log')
    fig.savefig(hp.plot_dir() + title + ".png", dpi = 300)
    fig.savefig(hp.plot_dir() + title + ".eps")

specs = {
        'nStSamp' : 5,
        'enGapSizes' : [6,8,10,12,14,16,18,20],
        'wVals' : [3.33, 3.42, 3.51, 3.6, 3.69, 3.78, 3.87]
        }

plot("GapCloseBBH3D", "locBBH3D_m1.1", specs, False, True)
#plot("LocBBH3D_intra1.1", "locBBH3D_intra1.1", specs)
