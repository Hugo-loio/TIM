import helper as hp
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

plt.style.use('science')

def shift_matrix_ticks(shape, ax):
    #x axis
    xticks = ax.get_xticks()[ax.get_xticks().tolist().index(0):]
    max_index = len(xticks)
    for i in range(1, len(xticks)):
        xticks[i] -= 1
        if(xticks[i] > shape[0]):
            max_index = i
    xticks = xticks[:max_index]
    xlabels = [str(int(tick + 1)) for tick in xticks]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)
    #y axis
    yticks = ax.get_yticks()[ax.get_yticks().tolist().index(0):]
    max_index = len(yticks)
    for i in range(1, len(yticks)):
        yticks[i] -= 1
        if(yticks[i] > shape[0]):
            max_index = i
    yticks = yticks[:max_index]
    ylabels = [str(int(tick + 1)) for tick in yticks]
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)

def ldos(fname, name, show: bool):
    data = hp.readfile(fname)
    x = data[0::3,0]
    y = data[1::3,0]
    zVec = np.average(data[2::3], axis = 1)
    #print(zVec)
    zMat = np.empty([int(x[-1])+1, int(y[-1])+1])
    zMin = 0
    zMax = 0.01
    for i in range(len(x)):
        z = zVec[int(y[i])*int(x[-1] + 1) + int(x[i])]
        zMat[int(x[i]),int(y[i])] = z
        if(z < zMin):
            zMin = z
        elif(z > zMax):
            zMax = z

    fig, ax = plt.subplots();
    im = ax.imshow(zMat, vmin = zMin, vmax = zMax, cmap = 'cividis')
    fig.colorbar(im , ax = ax)
    shift_matrix_ticks(zMat.shape, ax)

    ax.set(xlabel = r'$x$', ylabel = r'$y$')

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 300, bbox_inches='tight')
    if(show):
        plt.show()

def plots(fileNames, names, show):
    for i in range(len(fileNames)):
        ldos(fileNames[i] + ".dat", names[i], show)

weights = ['1', '2.6', '3.4']
fileNames = ['ldosSOTAI_L10_w' + w + '_E0_nMu1024_m1.1' for w in weights]
#plots(fileNames, fileNames, False)

weights = ['2.6', '3.4']
fileNames = ['ldosSOTAI_L10_w' + w + '_E0_nMu2048_m1.1' for w in weights]
#plots(fileNames, fileNames, False)

weights = ['2.6']
fileNames = ['ldosSOTAI_L10_w' + w + '_E0_nMu4096_m1.1' for w in weights]
#plots(fileNames, fileNames, False)

weights = ['1', '2.6', '3.4']
fileNames = ['ldosSOTAI_L40_w' + w + '_E0_nMu1024_m1.1' for w in weights]
#plots(fileNames, fileNames, False)

weights = ['1', '2.6', '3.4']
fileNames = ['ldosSOTAI_L60_w' + w + '_E0_nMu1024_m1.1' for w in weights]
#plots(fileNames, fileNames, False)

weights = ['5','8']
fileNames = ['ldosSOTAI_L60_w' + w + '_E0_nMu1024_m1.1' for w in weights]
plots(fileNames, fileNames, False)

weights = ['2.6']
fileNames = ['ldosSOTAI_L10_w' + w + '_E0_range0.01_m1.1_diag' for w in weights]
#plots(fileNames, fileNames, False)

weights = ['2.6']
fileNames = ['ldosSOTAI_L10_w' + w + '_E0_range0.05_m1.1_diag' for w in weights]
#plots(fileNames, fileNames, False)
