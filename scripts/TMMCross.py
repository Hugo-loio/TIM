import helper as hp
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy import optimize
import numpy as np

plt.style.use('science')


def fitFunc(x, a, b):
    return a*x + b

def chi2(x, y, yerr, a, b):
    return (1/(len(y) - 2))*sum(((y - fitFunc(x, a, b))**2)/yerr**2)

def min(num1, num2):
    if(num1 < num2):
        return num1
    else:
        return num2

def findCross(name1, name2):
    data1 = hp.readfile(name1 + ".dat")
    data2 = hp.readfile(name2 + ".dat")

    closeLambda = []
    closeW = []

    for i in range(min(len(data1[0]), len(data2[0]))):
        avgLambda = (data1[1][i] + data2[1][i])/2
        err = data1[2][i] + data2[2][i]
        if(abs(data1[1][i] - data2[1][i])/avgLambda < err):
            closeLambda.append(avgLambda)
            closeW.append(data1[0][i])

    '''
    if(len(closeLambda) == 0):
        for i in range(min(len(data1[0]), len(data2[0]))-1):
            avgLambda = (data1[1][i] + data2[1][i])/2
            avgLambda2 = (data1[1][i+1] + data2[1][i+1])/2
            diffLambda = (data1[1][i] - data2[1][i])/2
            diffLambda2 = (data1[1][i+1] - data2[1][i+1])/2
            if(diffLambda*diffLambda2 < 0):
                closeLambda.append(avgLambda)
                closeLambda.append(avgLambda2)
                closeW.append(data1[0][i])
                closeW.append(data1[0][i+1])
                '''

    #print(closeLambda)
    #print(np.average(closeLambda))
    #print(closeW)
    return [np.average(closeW), np.average(closeLambda), np.std(closeW)/np.sqrt(len(closeW)), np.std(closeLambda)/np.sqrt(len(closeLambda))]


def plotTMM(name, label, ax):
    data = hp.readfile(name + ".dat")
    data = data[:, data[0,:].argsort()]
    ax.plot(data[0], data[1], label = label, linestyle='-', linewidth = 0.5)


def plotConstL(name, fileNames, labels, crossVals, show = True):
    fig, ax = plt.subplots()

    for i in range(0, len(fileNames)):
        plotTMM(fileNames[i], labels[i], ax)

    ax.set(xlabel = r'$W$', ylabel = r'$\Lambda$')
    plt.yscale('log')
    ax.margins(x = 0)
    ax.errorbar(crossVals[0], crossVals[1], xerr = crossVals[2], yerr = crossVals[3], linewidth = 0.5,  capsize = 1, capthick = 0.5)
    #ax.legend(fontsize = 6, ncol = 2, bbox_to_anchor=(0.4,0.69))
    ax.legend(fontsize = 6, ncol = 2)

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 300)
    fig.savefig(hp.plot_dir() + name + ".eps")
    if(show):
        plt.show()

constLVals = ["4", "6", "8", "10", "12", "14"]
constLNames = ["tmmBBH3D_E0_L" + val + "_d2_m1.1_cross" for val in constLVals]
constLLabels = [r'$L_{x/y} = $ ' + val for val in constLVals]


crossValues = []
for i in range(len(constLVals) - 1):
    crossValues.append(findCross(constLNames[i], constLNames[i+1]))
crossValues = np.array(crossValues).T

plotConstL("tmmBBH3D_E0_d2_m1.1_cross", constLNames, constLLabels, crossValues, False)

fig, ax = plt.subplots()
lArray = np.array([int(val) for val in constLVals[:-1]])
x = 1/lArray
sort = np.argsort(x)
x = x[sort]
y = crossValues[0][sort]
yerr = crossValues[2][sort]
ax.errorbar(x, y, yerr = yerr, capsize = 1, linestyle="", linewidth = 0.5, capthick = 0.5)
ax.set(xlabel = r'$\frac{1}{L}$', ylabel = r'$W_{\times}$')
ax.set_xlim([0,0.3])
ax.set_ylim([3.5,4.05])
#plt.xscale('log')
#plt.yscale('log')
w0 = []
for i in range(3):
    x_temp = x[0:len(x)-i]
    y_temp = y[0:len(y)-i]
    yerr_temp = yerr[0:len(y)-i]
    popt, pcov = optimize.curve_fit(fitFunc, x_temp, y_temp, sigma = yerr_temp, absolute_sigma = True)
    perr = np.sqrt(np.diag(pcov))
    xFit = np.linspace(0, x_temp[-1], 100)
    ax.plot(xFit, fitFunc(xFit, popt[0], popt[1]), linewidth = 0.5)
    w0.append(popt[1])
#text = r'$ W_{\times}(0) = $'+ str(round(popt[1], 2)) + r'$ \pm $' + str(round(perr[1], 2))  + '\n'  r'$\chi^2_r = $' + str(round(chi2(x,y,yerr,popt[0], popt[1]), 2))
#ax.text(0.2, 3.65, text, ha = 'center')
w0_avg = np.average(w0)
w0_err = np.amax(w0-w0_avg)
text = r'$ W_{\times}(0) = $'+ str(round(w0_avg, 2)) + r'$ \pm $' + str(round(w0_err, 2))
ax.text(0.2, 3.65, text, ha = 'center')

#plt.show()
fig.savefig(hp.plot_dir() + "tmmBBH3D_E0_d2_m1.1_cross_L.png", dpi = 300)
fig.savefig(hp.plot_dir() + "tmmBBH3D_E0_d2_m1.1_cross_L.eps")

