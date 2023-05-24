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

def getEn(data):
    nsamp = int(len(data[1:])/nmax)
    res = np.empty((len(data[0]), nsamp, nmax))
    for i in range(len(data[0])):
        for j in range(nsamp):
            en = np.array(data[1+j*nmax:1+(j+1)*nmax,i])
            res[i,j] = en
    return res

def plotEnGap(name, fileName, show):
    data = hp.readfile(fileName + ".dat")
    w = data[0]
    nsamp = int(len(data[1:])/nmax)
    lsr = []
    lsrErr = []
    nums = np.arange(6,nmax+0.5,2)
    ls = getSpacings(data)
    aux = int((nmax-2)/2)
    ls1 = ls[:,:,:aux+3]
    ls2 = ls[:,:,aux+4:]
    ls = np.concatenate((ls1,ls2), axis = 2)
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

    avg_ratios = np.average(ratios, axis = 1)
    avg_spacings = np.average(ls[:,:,int((nmax-2)/2):], axis = 1)
    avg_max = np.average(max_ls, axis = 1)
    avg_min = np.average(min_ls, axis = 1)
    #print(np.sort(avg_spacings[0]))
    #print(avg_spacings[10])
    #print(avg_ratios[10])
    #print(avg_max[10])
    #print(avg_min[10])
    #en = getEn(data)
    #print(en[0,0,:])
    #print(en[2,0,:])
    #print(en[2,1,:])
    #print(en[2,2,:])
    #print(data[0][2])

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

    fig.savefig(hp.plot_dir() + name + "_v2.png", dpi = 300, bbox_inches='tight')
    if(show):
        plt.show()
        
    fig, ax = plt.subplots()
    ax.axhline(y = 0.53, color = 'black', linestyle = '--', linewidth = 0.5)
    ax.axhline(y = 0.386, color = 'black', linestyle = '--', linewidth = 0.5)
    ax.set(xlabel = r'$W$', ylabel = r'LSR')
    for n in [50, 30, 10]:
        i = nums.tolist().index(n)
        ax.errorbar(w[1:], lsr[i][1:], capsize = 2, linewidth = 0.5, capthick = 0.5, label = r'$n = ' + str(n) + '$')
    ymin, ymax = ax.get_ylim()
    ax.set_ylim([ymin, 1.01*ymax])
    #hp.totaiPhases(ax, 0.95)
    hp.totaiPhases(ax, 0.97)
    #fig.set_size_inches(1.7,1.3)
    fig.set_size_inches(4,3)
    #ax.legend(fontsize = 7, bbox_to_anchor = (0.41,0.55, 0.16, 0.1), alignment = 'center')
    ax.legend(bbox_to_anchor = (0.7,0.5), loc = 'center')
    ax.yaxis.labelpad = 0
    ax.tick_params(axis='y', which='major', pad=0.5)
    fig.savefig(hp.plot_dir() + name + ".png", dpi = 300, bbox_inches='tight')
    fig.savefig(hp.plot_dir() + name + ".pdf", bbox_inches='tight', pad_inches = 0.01)


sizes = ['20']
fileNames = ['spectrumE0BBH3D_L' + size + '_m1.1' for size in sizes]
for fileName in fileNames:
    plotEnGap(fileName + "LSR", fileName, False)
