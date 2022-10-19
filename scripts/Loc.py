import helper as hp
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

plt.style.use('science')

def plotCurve(ax, x, y, label, errors):
    if(errors):
        ax.errorbar(x, np.average(y, axis=0), label = label, yerr = np.std(y, axis = 0)/np.sqrt(len(y)) , capsize = 5, linestyle='-', linewidth = 0.5)
    else:
        ax.plot(x, np.average(y, axis=0), label = label, linestyle='-', linewidth = 0.5)

def fitFunc(x, a, b):
    return b*(1/x)**a

def plot(title, fileName, specs, errors = False, show = False):
    data = hp.readfile(fileName + ".dat")
    step = specs['nStSamp']*2 + 1

    if(len(specs['iprSizes']) > 0 and len(specs['iprNStates']) > 0):
        fig, ax = plt.subplots()
        for i in specs['iprSizes']:
            dataPlot = data[:,np.where(data[0] == i)[0]] 
            for e in specs['iprNStates']:
                start = int(2 + 2*(specs['nStSamp'] - e/10))
                label = r'$L = $ ' + str(i) + r' $ n = $ ' + str(e)
                plotCurve(ax, dataPlot[1], dataPlot[start::step], label, errors)

        ax.set(xlabel = r'$W$', ylabel = r'IPR')
        ax.legend(fontsize = 6)

        fig.savefig(hp.plot_dir() + title + "IPR.png", dpi = 300)
        fig.savefig(hp.plot_dir() + title + "IPR.eps")
        if(show):
            plt.show()
        plt.close()

    if(len(specs['lsrSizes']) > 0 and len(specs['lsrNStates']) > 0):
        fig, ax = plt.subplots()
        for i in specs['lsrSizes']:
            dataPlot = data[:,np.where(data[0] == i)[0]] 
            for e in specs['lsrNStates']:
                start = int(3 + 2*(specs['nStSamp'] - e/10))
                label = r'$L = $ ' + str(i) + r' $ n = $ ' + str(e)
                plotCurve(ax, dataPlot[1], dataPlot[start::step], label, errors)

        ax.set(xlabel = r'$W$', ylabel = r'LSR')
        ax.axhline(y = 0.53, color = 'black', linestyle = '--', linewidth = 0.5)
        ax.axhline(y = 0.386, color = 'black', linestyle = '--', linewidth = 0.5)
        ax.legend(fontsize = 6)

        fig.savefig(hp.plot_dir() + title + "LSR.png", dpi = 300)
        fig.savefig(hp.plot_dir() + title + "LSR.eps")
        if(show):
            plt.show()
        plt.close()

    if(len(specs['enGapSizes']) > 0):
        fig, ax = plt.subplots()
        for i in specs['enGapSizes']:
            dataPlot = data[:,np.where(data[0] == i)[0]] 
            start = int(2 + 2*specs['nStSamp'])
            label = r'$L = $' + str(i)
            plotCurve(ax, dataPlot[1], dataPlot[start::step], label, errors)

        ax.set(xlabel = r'$W$', ylabel = r'Gap')
        ax.legend(fontsize = 6)

        fig.savefig(hp.plot_dir() + title + "Gap.png", dpi = 300)
        fig.savefig(hp.plot_dir() + title + "Gap.eps")
        if(show):
            plt.show()
        plt.close()

    if(len(specs['fractalSizes']) > 0 and len(specs['fractalNStates']) > 0):
        fig, ax = plt.subplots()
        dataPlot = data[:,np.where(data[0] == specs['fractalSizes'][0])[0]] 
        for i in range(1, len(specs['fractalSizes'])):
            dataPlot = np.concatenate((dataPlot, data[:,np.where(data[0] == specs['fractalSizes'][i])[0]]), axis=1)

        nPoints = len(np.where(dataPlot[0] == dataPlot[0][0])[0])
        for n in specs['fractalNStates']:
            start = int(2 + 2*(specs['nStSamp'] - n/10))
            vals = np.average(dataPlot[start::step], axis = 0)
            errs = np.std(dataPlot[start::step], axis = 0)/np.sqrt(len(dataPlot[start::step,0]))
            d = np.zeros(nPoints)
            dErr = np.zeros(nPoints)
            for i in range(0, nPoints):
                xFit = dataPlot[0][i::nPoints]
                yFit = vals[i::nPoints]
                yFitErr = errs[i::nPoints]
                if(0 in yFitErr):
                    popt, pcov = optimize.curve_fit(fitFunc, xFit, yFit)
                else:
                    popt, pcov = optimize.curve_fit(fitFunc, xFit, yFit, sigma = yFitErr)
                d[i] = popt[0]
                dErr[i] = np.sqrt(np.diag(pcov))[0]

            ax.plot(dataPlot[1][0:nPoints], d, linestyle = '-', label = "n = " + str(n))

        ax.set(xlabel = r'$W$', ylabel = r'$D_2$')
        ax.legend(fontsize = 6)
        fig.savefig(hp.plot_dir() + title + "Fractal.png", dpi = 300)
        fig.savefig(hp.plot_dir() + title + "Fractal.eps")
        if(show):
            plt.show()
        plt.close()

specs = {
        'nStSamp' : 5,
        'iprSizes' : [260],
        'iprNStates' : [50,30,10],
        'fractalNStates' : [50,30,10],
        'fractalSizes' : [100,120,140,160,180,200,220,260],
        'lsrSizes' : [260],
        'lsrNStates' : [50,30,10],
        'enGapSizes' : [260]
        }

plot("LocSOTAI_m1.1", "locSOTAI_m1.1", specs)

specs = {
        'nStSamp' : 5,
        'iprSizes' : [12, 14, 16],
        'iprNStates' : [50,10],
        'fractalNStates' : [50,40,30,20,10],
        'fractalSizes' : [4,6,8,10,12,14,16],
        'lsrSizes' : [12, 14, 16],
        'lsrNStates' : [50,10],
        'enGapSizes' : [12, 14, 16]
        }

#plot("LocBBH3D_m1.1", "locBBH3D_m1.1", specs)
#plot("LocBBH3D_intra1.1", "locBBH3D_intra1.1", specs)
