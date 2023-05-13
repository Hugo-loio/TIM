import helper as hp
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import gridspec
import matplotlib.colors as cls
import numpy as np

plt.style.use('science')

nmax = 50

def getSpacings(data):
    nsamp = int(len(data[1:])/nmax)
    aux = int(nmax/2-1)
    res = np.empty((len(data[0]), nsamp, nmax-2))
    for i in range(len(data[0])):
        for j in range(nsamp):
            en1 = np.array(data[1+j*nmax:(j+1)*nmax,i])
            en2 = np.array(data[2+j*nmax:1+(j+1)*nmax,i])
            ls = np.empty(nmax-2)
            ls[:aux] = en2[:aux] - en1[:aux]
            ls[aux:] = en2[aux+1:] - en1[aux+1:]
            res[i,j] = ls
    return res

def plotEnGap(name, fileName, show):
    data = hp.readfile(fileName + ".dat")
    w = data[0]
    nsamp = int(len(data[1:])/nmax)
    lsr = []
    lsrErr = []
    nums = np.arange(6,nmax+0.5,2)
    ls = getSpacings(data)
    max_ls = np.maximum(ls[:,:,int((nmax-2)/2):-1], ls[:,:,int((nmax-2)/2+1):])
    min_ls = np.minimum(ls[:,:,int((nmax-2)/2):-1], ls[:,:,int((nmax-2)/2+1):])
    max_ls[max_ls == 0] = 1
    min_ls[min_ls == 0] = 1
    ratios = min_ls/max_ls
    for n in nums:
        n_eff = int(n/2)
        lsrs = np.average(ratios[:,:,:n_eff-2], axis = 2)
        lsr.append(np.average(lsrs, axis = 1))
        lsrErr.append(np.std(lsrs, axis = 1)/np.sqrt(nsamp))



    spec = gridspec.GridSpec(ncols=2, nrows=1, width_ratios=[40, 1])
    fig = plt.figure()
    ax = fig.add_subplot(spec[0])
    ax.set(xlabel = r'$W$', ylabel = r'LSR')

    black = cls.to_rgba_array('black')[0]
    color = np.array([float(num) for num in cm.tab10(0)])
    black = cls.to_rgba_array('black')[0]
    white = cls.to_rgba_array('white')[0]
    colors = [color*0.25 + 0.75*white, color,color*0.5 + 0.50*black]
    nodes = [0.0, 0.5, 1]
    cmap = cls.LinearSegmentedColormap.from_list("mycmap", list(zip(nodes,colors)))
    vmin = np.amin(nums)
    vmax = np.amax(nums)
    ax_cmap = fig.add_subplot(spec[1])
    pos = ax_cmap.get_position()
    ax_cmap.set_position([0.95*pos.x0, pos.y0, pos.width, pos.height])
    norm = cls.Normalize(vmin = vmin, vmax = vmax)
    cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax = ax_cmap, orientation = 'vertical')
    #cbar.ax.set_title(kwargs['cbar_title'])
    cbar.set_label('n')

    #plt.yscale('log')
    #plt.xscale('log')
    #ax.margins(x = 0)

    for i,n in enumerate(nums[0:]):
        color = cmap((n-vmin)/(vmax-vmin))
        ax.errorbar(w[1:], lsr[i][1:], color = color, capsize = 2, linewidth = 0.5, capthick = 0.5)
    ax.axhline(y = 0.53, color = 'black', linestyle = '--', linewidth = 0.5)
    ax.axhline(y = 0.386, color = 'black', linestyle = '--', linewidth = 0.5)
    hp.totaiPhases(ax,0.95)
    #ax.legend(fontsize = 4, ncol = 2, title = "n")

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 300, bbox_inches='tight')
    if(show):
        plt.show()

sizes = ['20']
fileNames = ['spectrumE0BBH3D_L' + size + '_m1.1' for size in sizes]
for fileName in fileNames:
    plotEnGap(fileName + "LSR", fileName, False)
