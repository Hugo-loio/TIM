import helper as hp
import matplotlib.pyplot as plt
import numpy as np

plt.style.use('science')

def plotCurve(ax, x, y, label, errors):
    if(errors):
        ax.errorbar(x, np.average(y, axis=0), label = label, yerr = np.std(y, axis = 0)/np.sqrt(len(y)) , capsize = 5, linestyle='-', linewidth = 0.5)
    else:
        ax.plot(x, np.average(y, axis=0), label = label, linestyle='-', linewidth = 0.5)

def plot(title, fileName, specs, errors = False, show = False):
    data = hp.readfile(fileName + ".dat")
    step = specs['nStSamp']*2 + 1

    if(len(specs['iprSizes']) > 0 and len(specs['iprNStates']) > 0):
        fig, ax = plt.subplots()
        for i in specs['iprSizes']:
            dataPlot = data[:,np.where(data[0] == i)[0]] 
            for e in specs['iprNStates']:
                start = int(2 + 2*(specs['nStSamp'] - e/10))
                label = r'$L = $' + str(i) + r' $ n = $' + str(e)
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
                label = r'$L = $' + str(i) + r' $ n = $' + str(e)
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

specs = {
        'nStSamp' : 5,
        'iprSizes' : [40,60],
        'iprNStates' : [50,30,10],
        'fractalNStates' : [50,40,30,20,10],
        'fractalSizes' : [50,60],
        'lsrSizes' : [40,60],
        'lsrNStates' : [50,30,10],
        'enGapSizes' : [40,60]
        }

plot("LocSOTAI_m1.1", "locSOTAI_m1.1", specs)

specs = {
        'nStSamp' : 5,
        'iprSizes' : [4, 6, 8],
        'iprNStates' : [50,30,10],
        'fractalNStates' : [50,40,30,20,10],
        'fractalSizes' : [4, 6, 8],
        'lsrSizes' : [4, 6, 8],
        'lsrNStates' : [50,30,10],
        'enGapSizes' : [4, 6, 8]
        }

plot("LocBBH3D_m1.1", "locBBH3D_m1.1", specs)
#plot("LocBBH3D_intra1.1", "locBBH3D_intra1.1", specs)
