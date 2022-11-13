import helper as hp
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

#plt.rcParams.update({'font.size': 14})
plt.style.use('science')

def fitFunc(x, a, b):
    return a*x + b

def plot(name):
    fig, ax = plt.subplots()
    data = hp.readfile(name + ".dat")

    data = data[:, data[0,:].argsort()]
    qyz = data[1::3]
    qxz = data[2::3]
    qxy = data[3::3]
    q = 8*np.absolute(np.multiply(np.multiply(qyz,qxz),qxy))

    x = 1/data[0]
    sort = np.argsort(x)
    x = x[sort]
    qavg = np.average(q, axis = 0)
    qavg = qavg[sort]
    qerr = np.std(q, axis = 0)/np.sqrt(len(q))
    qerr = qerr[sort]

    #ax.errorbar(x, qavg, yerr = qerr, capsize = 4, linestyle='none', linewidth = 0.5, capthick = 0.5, marker = '_')
    ax.errorbar(x, qavg, yerr = qerr, capsize = 4, linestyle='none', marker = '_')

    ax.set(xlabel = r'$\frac{1}{L}$', ylabel = r'$Q$')

    nFitPoints = [5,10,len(x)]
    bFit = []
    for fitEnd in nFitPoints:
        xfit = x[:fitEnd]
        yfit = qavg[:fitEnd]
        yfiterr = qerr[:fitEnd] 

        popt, pcov = optimize.curve_fit(fitFunc, xfit, yfit, sigma = yfiterr, absolute_sigma = True)
        xdata = np.linspace(0, x[fitEnd-1], 100)
        ax.plot(xdata, fitFunc(xdata, popt[0], popt[1]), linewidth = 0.5)
        bFit.append(popt[1])
        
    #print(np.average(bFit))
    #print(bFit)

    '''
    fitStart = 0
    xfit = x[fitStart:]
    yfit = qavg[fitStart:]
    yfiterr = qerr[fitStart:] 

    popt, pcov = optimize.curve_fit(fitFunc, xfit, yfit, sigma = yfiterr, absolute_sigma = True)
    perr = np.sqrt(np.diag(pcov))
    params = r'$ Q(L) = \frac{a}{L} + b$' '\n'  r'$a = $' + str(round(popt[0], 2)) + r'$ \pm $' + str(round(perr[0],2)) +  '\n' + r'$b = $' + str(round(popt[1], 2)) + r'$ \pm $' + str(round(perr[1],2))

    xdata = np.linspace(0, x[-1], 100)
    ax.plot(xdata, fitFunc(xdata, popt[0], popt[1]), label = params, color = 'orange')
    '''

    #plt.legend(loc= 'upper right')
    text = r'$ Q(0) = $'+ str(round(np.average(bFit), 2)) + r'$ \pm $' + str(round(np.std(bFit)/np.sqrt(len(bFit)),2))
    xmax = ax.get_xlim()[1]
    ymax = ax.get_ylim()[1]
    ax.text(0.7*xmax, 0.9*ymax, text, ha = 'center')
    fig.savefig(hp.plot_dir() + name + ".png", dpi = 300)
    fig.savefig(hp.plot_dir() + name + ".eps")
    #plt.show()
    #plt.close()

plot("constantDisorderBBH3Dquad_intra1.1_w2.8")
plot("constantDisorderBBH3Dquad_intra1.1_w3")
plot("constantDisorderBBH3Dquad_intra1.1_w3.2")
plot("constantDisorderBBH3Dquad_intra1.1_w3.4")
plot("constantDisorderBBH3Dquad_intra1.1_w3.6")
plot("constantDisorderBBH3Dquad_intra1.1_w4")
plot("constantDisorderBBH3Dquad_intra1.1_w9")
