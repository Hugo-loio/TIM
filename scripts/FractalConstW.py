import helper as hp
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

plt.style.use('science')

def fitFunc(x, a, b):
    return b*(1/x)**a

def plot(name, fileNames, labels, sizes, show = True):
    fig, ax = plt.subplots()
    for e in range(len(fileNames)):
        data = hp.readfile(fileNames[e] + ".dat")
        #data = data[:, data[0,:].argsort()]

        if(len(sizes) > 0):
            dataPlot = data[:,np.where(data[0] == sizes[0])[0]]
            for i in range(1, len(sizes)):
                dataPlot = np.concatenate((dataPlot, data[:,np.where(data[0] == sizes[i])[0]]), axis=1)

            nPoints = len(np.where(dataPlot[0] == dataPlot[0][0])[0])
            vals = np.average(dataPlot[2:], axis = 0)
            errs = np.std(dataPlot[2:], axis = 0)/np.sqrt(len(dataPlot[2:,0]))

            d = np.zeros(nPoints)
            dErr = np.zeros(nPoints)

            for i in range(0, nPoints):
                xFit = dataPlot[0][i::nPoints]
                yFit = vals[i::nPoints]
                yFitErr = errs[i::nPoints]
                popt, pcov = optimize.curve_fit(fitFunc, xFit, yFit, sigma = yFitErr)

                d[i] = popt[0]
                dErr[i] = np.sqrt(np.diag(pcov))[0]

            #ax.errorbar(dataPlot[1][0:nPoints], d, yerr = dErr, capsize = 5, linestyle = '-', label = labels[e])
            ax.plot(dataPlot[1][0:nPoints], d, linestyle = '-', label = labels[e])

    ax.set(xlabel = r'$E$', ylabel = r'$D_2$')
    ax.legend()

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 300)
    fig.savefig(hp.plot_dir() + name + ".eps")

    if(show):
        plt.show()
    plt.close()


weights = ["3.2", "4.6"]
sizes = [120, 140, 160, 164, 168, 172, 180, 184, 188, 192, 196, 200, 204, 208, 212, 216, 220, 224, 228, 232, 236, 240]
names = ["fractalSOTAI_w" + w + "_m1.1" for w in weights]
labels = ["W = " + w for w in weights]
#plot("FractalConstWSOTAI", names, labels, sizes, False)

weights = ["3.4", "5"]
sizes = [4, 6, 8, 10, 12, 14, 16]
names = ["fractalBBH3D_w" + w + "_m1.1" for w in weights]
labels = ["W = " + w for w in weights]
plot("FractalConstWBBH3D", names, labels, sizes, False)

