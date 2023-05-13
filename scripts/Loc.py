import helper as hp
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

plt.style.use('science')

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Times New Roman",
    "font.size": 10
})

def plotCurve(ax, x, y, label, errors):
    if(errors):
        ax.errorbar(x, np.average(y, axis=0), label = label, yerr = np.std(y, axis = 0)/np.sqrt(len(y)) , capsize = 5, linestyle='-', linewidth = 0.5)
    else:
        ax.plot(x, np.average(y, axis=0), label = label, linestyle='-')

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
                if(len(specs['iprSizes']) == 1):
                    label = r'$ n = $ ' + str(e)
                else:
                    label = r'$L = $ ' + str(i) + r' $ n = $ ' + str(e)
                plotCurve(ax, dataPlot[1], dataPlot[start::step], label, errors)

        ax.set(xlabel = r'$W$', ylabel = r'IPR')
        ax.margins(x = 0)
        if(specs['phases'] == 0):
            hp.sotaiPhases(ax)
            #ax.legend(fontsize = 6)
            ax.legend(fontsize = 6, ncol = 1, bbox_to_anchor=(0.4,0.7))
        else:
            hp.totaiPhases(ax, 0.85)
            ax.legend(fontsize = 6, bbox_to_anchor = (0.45,0.9))

        ax.tick_params(axis='y', which='major', pad=0.5)
        ax.yaxis.labelpad = 0
        fig.set_size_inches(1.7,1.3)
        fig.savefig(hp.plot_dir() + title + "IPR.png", bbox_inches = 'tight', dpi = 300)
        fig.savefig(hp.plot_dir() + title + "IPR.pdf", bbox_inches = 'tight', pad_inches = 0)
        if(show):
            plt.show()
        plt.close()

    if(len(specs['lsrSizes']) > 0 and len(specs['lsrNStates']) > 0):
        fig, ax = plt.subplots()
        for i in specs['lsrSizes']:
            dataPlot = data[:,np.where(data[0] == i)[0]] 
            for e in specs['lsrNStates']:
                start = int(3 + 2*(specs['nStSamp'] - e/10))
                if(len(specs['lsrSizes']) == 1):
                    label = r'$ n = $ ' + str(e)
                else:
                    label = r'$L = $ ' + str(i) + r' $ n = $ ' + str(e)
                plotCurve(ax, dataPlot[1], dataPlot[start::step], label, errors)

        ax.set(xlabel = r'$W$', ylabel = r'LSR')
        ax.axhline(y = 0.53, color = 'black', linestyle = '--', linewidth = 0.5)
        ax.axhline(y = 0.386, color = 'black', linestyle = '--', linewidth = 0.5)
        ax.margins(x = 0)
        if(specs['phases'] == 0):
            hp.sotaiPhases(ax, 0.95)
        else:
            ymin, ymax = ax.get_ylim()
            ax.set_ylim([ymin, 1.03*ymax])
            hp.totaiPhases(ax, 0.95)
            ax.legend(fontsize = 7, bbox_to_anchor = (0.4,0.75))
        #ax.legend(fontsize = 7)
        ax.tick_params(axis='y', which='major', pad=0.5)
        ax.yaxis.labelpad = 0

        fig.set_size_inches(1.7,1.3)
        fig.savefig(hp.plot_dir() + title + "LSR.png", bbox_inches = 'tight', dpi = 300)
        fig.savefig(hp.plot_dir() + title + "LSR.pdf", bbox_inches = 'tight', pad_inches = 0)
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
        ax.margins(x = 0)
        if(len(specs['enGapSizes']) > 1): 
            ax.legend(fontsize = 6)
        if(specs['phases'] == 0):
            hp.sotaiPhases(ax)

        fig.set_size_inches(1.7, 1.3)
        hp.totaiPhases(ax, 0.8)
        fig.savefig(hp.plot_dir() + title + "Gap.png", bbox_inches = 'tight', dpi = 300, pad_inches = 0.01)
        fig.savefig(hp.plot_dir() + title + "Gap.pdf", bbox_inches = 'tight', pad_inches = 0.01)
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
        ax.axhline(y = 3, color = 'black', linestyle = '--', linewidth = 0.5)
        ax.margins(x = 0)
        if(specs['phases'] == 0):
            hp.sotaiPhases(ax)
            ax.legend(fontsize = 6, bbox_to_anchor=(0.71,0.8))
        else:
            #ax.legend(fontsize = 6)
            ymin, ymax = ax.get_ylim()
            ax.set_ylim([ymin, 1.1*ymax])
            hp.totaiPhases(ax, 0.85)
            ax.legend(fontsize = 7, bbox_to_anchor = (0.4,0.60))

        ax.tick_params(axis='y', which='major', pad=0.5)
        ax.yaxis.labelpad = 0

        fig.set_size_inches(1.7,1.3)
        fig.savefig(hp.plot_dir() + title + "Fractal.png", bbox_inches = 'tight', dpi = 300)
        fig.savefig(hp.plot_dir() + title + "Fractal.pdf", bbox_inches = 'tight', pad_inches = 0)
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
        'enGapSizes' : [260],
        'phases' : 0
        }

#plot("LocSOTAI_m1.1", "locSOTAI_m1.1", specs)

specs = {
        'nStSamp' : 5,
        'iprSizes' : [20],
        'iprNStates' : [50, 30, 10],
        'fractalNStates' : [50,30,10],
        'fractalSizes' : [10,12,14,16,18,20],
        'lsrSizes' : [20],
        'lsrNStates' : [50, 30,10],
        'enGapSizes' : [20],
        'phases' : 1
        }

plot("LocBBH3D_m1.1", "locBBH3D_m1.1", specs, False, False)
#plot("LocBBH3D_intra1.1", "locBBH3D_intra1.1", specs)
