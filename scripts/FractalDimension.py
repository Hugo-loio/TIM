import helper as hp
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

plt.style.use('science')

def fitFunc(x, a, b):
    return b*(1/x)**a

def plot(name, fileName, show = True):
    data = hp.readfile(fileName + ".dat")
    #data = data[:, data[0,:].argsort()]

    nPoints = len(np.where(data[0] == data[0][0])[0])
    vals = np.average(data[2:], axis = 0)
    errs = np.std(data[2:], axis = 0)/np.sqrt(len(data[2:,0]))

    d = np.zeros(nPoints)
    dErr = np.zeros(nPoints)

    for i in range(0, nPoints):
        xFit = data[0][i::nPoints]
        yFit = vals[i::nPoints]
        yFitErr = errs[i::nPoints]
        popt, pcov = optimize.curve_fit(fitFunc, xFit, yFit, sigma = yFitErr)

        d[i] = popt[0]
        dErr[i] = np.sqrt(np.diag(pcov))[0]

    fig, ax = plt.subplots()

    #ax.errorbar(data[1][0:nPoints], d, yerr = dErr, capsize = 5, linestyle = '-')
    ax.plot(data[1][0:nPoints], d, linestyle = '-')

    ax.set(xlabel = r'$W$', ylabel = r'$D_2$')
    #ax.legend()

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 300)

    if(show):
        plt.show()
    plt.close()


name = "iprSOTAI_n10_m1.1"
plot("FractalDimensionSOTAI_n10", name, False)

name = "iprSOTAI_n50_m1.1"
plot("FractalDimensionSOTAI_n50", name, False)
