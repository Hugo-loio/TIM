import helper as hp
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

def fitFunc(x, a, b):
    return a*(1/x) + b

def plot(name):
    fig, ax = plt.subplots()
    data = hp.readfile(name + ".dat")

    data = data[:, data[0,:].argsort()]
    qyz = data[1::3]
    qxz = data[2::3]
    qxy = data[3::3]
    q = 8*np.absolute(np.multiply(np.multiply(qyz,qxz),qxy))
    #print(q[0:,-1])
    print(np.sqrt(len(q)))

    ax.errorbar(data[0], np.average(q,axis=0), yerr = np.std(q, axis = 0)/np.sqrt(len(q)), capsize = 5, linestyle='-')

    ax.set(xlabel = r'$L$', ylabel = r'$Q$')

    fitStart = 5
    xfit = data[0][fitStart:]
    yfit = np.average(q, axis=0)[fitStart:]
    yfiterr = np.std(q, axis = 0)[fitStart:]/np.sqrt(len(q) - fitStart)

    popt, pcov = optimize.curve_fit(fitFunc, xfit, yfit, sigma = yfiterr)
    perr = np.sqrt(np.diag(pcov))
    params = r'$ f(x) = \frac{a}{x} + b \ , \  a = $' + str(round(popt[0], 2)) + r'$ \pm $' + str(round(perr[0],2)) +  r'$ \ , \ b = $' + str(round(popt[1], 2)) + r'$ \pm $' + str(round(perr[1],2))
    print(params)

    xdata = np.linspace(data[0][fitStart], data[0][-1], 100)
    ax.plot(xdata, fitFunc(xdata, popt[0], popt[1]), label = params)
    plt.legend()

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 200)
    plt.show()
    plt.close()

plot("constantDisorderBBH3Dquad_intra1.1_w3")
