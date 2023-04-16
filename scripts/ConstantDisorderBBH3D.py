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
    ax.errorbar(x, qavg, yerr = qerr, capsize = 1.5, linestyle='none', marker = '_', linewidth = 0.5, capthick = 0.5, markersize = 2.5, markeredgewidth = 0.5)

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

    #plt.legend(loc= 'upper right')
    q0 = np.average(bFit)
    q0Err = (np.absolute(bFit-q0)).max()
    text = r'$ Q(0) = $'+ str(round(q0, 2)) + r'$ \pm $' + str(round(q0Err,2))
    xmax = ax.get_xlim()[1]
    ymax = ax.get_ylim()[1]
    ax.text(0.55*xmax, 0.85*ymax, text, ha = 'center')

    fig.set_size_inches(1.7, 1.3)
    fig.savefig(hp.plot_dir() + name + ".png", dpi = 300, bbox_inches = 'tight', pad_inches = 0.01)
    fig.savefig(hp.plot_dir() + name + ".pdf", bbox_inches = 'tight', pad_inches = 0.01)
    return q0, q0Err
    #plt.show()
    #plt.close()

wVec = [2.8, 3, 3.2, 3.4, 3.6, 4, 9] 
q0Vec = []
q0ErrVec = []
for w in wVec:
    q0, q0Err = plot("constantDisorderBBH3Dquad_intra1.1_w" + str(w))
    q0Vec.append(q0)
    q0ErrVec.append(q0Err)

res = [wVec, q0Vec, q0ErrVec]
