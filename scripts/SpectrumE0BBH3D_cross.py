import helper as hp
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import gridspec
import matplotlib.colors as cls
from scipy import optimize
import numpy as np

plt.style.use('science')

def fitFunc(x, a, b):
    return a*x + b

def chi2(x, y, yerr, a, b):
    return (1/(len(y) - 2))*sum(((y - fitFunc(x, a, b))**2)/yerr**2)

def getEnGap(data):
    gap = []
    gapErr = []
    for i in range(len(data[0])):
        gaps = [data[j+1,i] - data[j,i] for j in range(1, len(data) -1) if (data[j+1,i] > 0 and data[j,i] < 0) ]
        gap.append(np.average(gaps))
        gapErr.append(np.std(gaps)/np.sqrt(len(gaps)))
    return np.array(gap), np.array(gapErr)

def getMLS(data, nums, nmax):
    mls = []
    mlsErr = []
    nsamp = len(data)/nmax
    midl = nmax/2-1
    for i in range(len(data[0])):
        ls = [data[j+1,i] - data[j,i] for j in range(1, len(data) -1) if (data[j+1,i]*data[j,i] > 0) ]
        mls.append([])
        mlsErr.append([])
        for n in nums:
            n2 =n/2-1
            lsN = np.array([ls[int(midl*i-n2):int(midl*i+n2)] for i in range(1, int(nsamp+1))]).flatten()
            mls[i].append(np.average(lsN))
            mlsErr[i].append(np.std(lsN)/np.sqrt(len(lsN)))
    return np.array(mls), np.array(mlsErr)

def getCross(fileNames, nums):
    x, xerr, y, yerr = [], [], [], []
    fig, ax = plt.subplots()
    ax.set(xlabel = r'$W$', ylabel = r'Energy')
    gaps, gapErrs, mlss, mlsErrs = [], [], [], []
    for fileName in fileNames:
        data = hp.readfile(fileName + ".dat")
        w = data[0]
        gap, gapErr = getEnGap(data)
        gaps.append(gap)
        gapErrs.append(gapErr)
        mls, mlsErr = getMLS(data, nums, nums[-1])
        mlss.append(mls)
        mlsErrs.append(mlsErr)

    for i,gap in enumerate(gaps):
        xs = []
        ys = []
        for e in range(len(gap)):
            if(abs(gap[e]-mlss[i][e,-1]) < gapErrs[i][e]):
                xs.append(w[e])
                ys.append(gap[e])
        x.append(np.average(xs))
        xerr.append(np.std(xs)/np.sqrt(len(xs)))
        y.append(np.average(ys))
        yerr.append(np.std(ys)/np.sqrt(len(ys)))

    spec = gridspec.GridSpec(ncols=2, nrows=2, width_ratios=[40, 1])
    fig = plt.figure()
    ax = fig.add_subplot(spec[:,0])
    ax.set(xlabel = r'$W$', ylabel = r'Energy')

    black = cls.to_rgba_array('black')[0]
    color1 = np.array([float(num) for num in cm.tab10(0)])
    color2 = np.array([float(num) for num in cm.tab10(3)])
    black = cls.to_rgba_array('black')[0]
    white = cls.to_rgba_array('white')[0]
    colors1 = [color1*0.25 + 0.75*white, color1,color1*0.5 + 0.50*black]
    colors2 = [color2*0.25 + 0.75*white, color2,color2*0.5 + 0.50*black]
    nodes = [0.0, 0.5, 1]
    cmap1 = cls.LinearSegmentedColormap.from_list("mycmap1", list(zip(nodes,colors1)))
    cmap2 = cls.LinearSegmentedColormap.from_list("mycmap2", list(zip(nodes,colors2)))
    L = np.array([int(size) for size in sizes])
    vmin = np.amin(L)
    vmax = np.amax(L)

    ax_cmap1 = fig.add_subplot(spec[0,1])
    pos = ax_cmap1.get_position()
    ax_cmap1.set_position([0.95*pos.x0, pos.y0, pos.width, pos.height])
    norm = cls.Normalize(vmin = vmin, vmax = vmax)
    cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap1), cax = ax_cmap1, orientation = 'vertical')
    cbar.set_label(r'$L$')
    cbar.ax.set_title('Gap', loc = 'left')

    ax_cmap2 = fig.add_subplot(spec[1,1])
    pos = ax_cmap2.get_position()
    ax_cmap2.set_position([0.95*pos.x0, 0.30*pos.y0, pos.width, pos.height])
    norm = cls.Normalize(vmin = vmin, vmax = vmax)
    cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap2), cax = ax_cmap2, orientation = 'vertical')
    cbar.set_label(r'$L$')
    cbar.ax.set_title('MLS', loc = 'left')
    ax.set_yscale('log')

    #plt.xscale('log')
    #ax.margins(x = 0)

    for i,gap in enumerate(gaps):
        color = cmap1((L[i]-vmin)/(vmax-vmin))
        ax.errorbar(w, gap, color = color, linewidth = 0.5)
        #ax.errorbar(w, gap, yerr = gapErrs[i], capsize = 2, linewidth = 0.5, capthick = 0.5, label = "Gap")

    for i,mls in enumerate(mlss):
        for e,num in enumerate(nums[-1:]):
            color = cmap2((L[i]-vmin)/(vmax-vmin))
            ax.errorbar(w, mls[:,e], color = color, linewidth = 0.5)
            #ax.errorbar(w, mls[:,e], yerr = mlsErr, capsize = 2, linewidth = 0.5, capthick = 0.5, label = "MLS")

    ax.errorbar(x, y, xerr = xerr, yerr = yerr, capsize = 2, linewidth = 0.5, capthick = 0.5, color = 'black', ls = '')
    #ax.legend(fontsize = 10, ncol = 1)
    #ax.legend()

    ax.yaxis.labelpad = 0
    fig.set_size_inches(3.3,2.5)
    ax.tick_params(axis='y', which='major', pad=0.5)
    fig.savefig(hp.plot_dir() + "spectrumE0BBH3D_m1.1_cross.png", dpi = 300, bbox_inches='tight', pad_inches = 0.01)
    fig.savefig(hp.plot_dir() + "spectrumE0BBH3D_m1.1_cross.pdf", bbox_inches='tight', pad_inches = 0.01)

    return np.array(x), np.array(xerr), np.array(y), np.array(yerr)

sizes = ['10', '12', '14', '16', '18', '20']
#sizes = sizes[0:1]
nums = [50]
fileNames = ['spectrumE0BBH3D_L' + size + '_m1.1_cross' for size in sizes]

x, xerr, y, yerr = getCross(fileNames, nums)

fig, ax = plt.subplots()
L = 1/np.array([int(size) for size in sizes])
sort = np.argsort(L)
L = L[sort]
x = x[sort]
xerr = xerr[sort]
ax.errorbar(L, x, yerr = xerr, capsize = 2, capthick = 0.5, linestyle="", linewidth = 0.5)
ax.set(xlabel = r'$\frac{1}{L}$', ylabel = r'$W_{\times}$')

w0 = []
w0_err = 0
for i in range(3):
    x_temp = L[0:len(x)-i]
    y_temp = x[0:len(y)-i]
    yerr_temp = xerr[0:len(y)-i]
    popt, pcov = optimize.curve_fit(fitFunc, x_temp, y_temp, sigma = yerr_temp, absolute_sigma = True)
    perr = np.sqrt(np.diag(pcov))
    xFit = np.linspace(0, x_temp[-1], 100)
    ax.plot(xFit, fitFunc(xFit, popt[0], popt[1]), linewidth = 0.5)
    w0.append(popt[1])
    w0_err = np.amax([perr[1], w0_err])

w0_avg = np.average(w0)
w0_err = np.amax([w0_err,np.amax(abs(w0-w0_avg))])
text = r'$ W_{\times}(0) = $'+ str(round(w0_avg, 2)) + r'$ \pm $' + str(round(1.5*w0_err, 2))
ax.text(0.04, 3.8, text, ha = 'center')

fig.set_size_inches(3.3,2.5)
fig.savefig(hp.plot_dir() + "spectrumE0BBH3D_m1.1_cross_fit.png", dpi = 300, bbox_inches='tight', pad_inches = 0.01)
fig.savefig(hp.plot_dir() + "spectrumE0BBH3D_m1.1_cross_fit.pdf", bbox_inches='tight', pad_inches = 0.01)
