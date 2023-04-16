import helper as hp
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

plt.style.use('science')

tran_num = np.array([[1.1, 2.55]])

def plot(name):
    fig, ax = plt.subplots()

    path_name = os.path.dirname(os.path.realpath(__file__)) + "/data/" + name + ".npy"
    with open(path_name, 'rb') as f:
        ren_gammas = np.load(f)

    path_name = os.path.dirname(os.path.realpath(__file__)) + "/data/" + name + "_params.npy"
    with open(path_name, 'rb') as f:
        axes = np.load(f)

    x = axes[0]
    y = axes[1]
    mat = np.amax(ren_gammas, axis = 2) 

    im = ax.pcolormesh(x, y, mat, shading = 'nearest', cmap = 'viridis', linewidth = 0, rasterized = True)
    cbar = fig.colorbar(im)
    cbar.set_label(r'$\max_{i}\{\gamma^\prime_i\}$')
    ax.set(xlabel = r'$W$', ylabel = r'$\gamma$')

    gamma_crit, w_crit = [], []
    for i,gamma in enumerate(y):
        for e,W in enumerate(x):
            if(e != 0):
                if(mat[i][e-1] >= 1 and mat[i][e] < 1):
                    w_crit.append((x[e-1] + x[e])/2)
                    gamma_crit.append(gamma)
                    break
            if(e == len(x)-1):
                print("Transition not found for gamma", gamma)

    ax.plot(w_crit,gamma_crit,'--', color = 'red', label = 'SCBA')
    ax.plot(tran_num[:,1], tran_num[:,0], 'o', label = r'$Q$', markersize = 6, color = 'blue')
    ax.legend(fontsize = 9, title = 'Transition')

    '''
    check_2d_ax_args(ax, kwargs)
    if 'cbar_title' in kwargs:
        cbar.ax.set_title(kwargs['cbar_title'])
    if 'cbar_label' in kwargs:
        cbar.set_label(kwargs['cbar_label'])
    ymin, ymax = ax.get_ylim()
    ax.set_ylim([ymin, ymax*1.2])
    ax.yaxis.labelpad = 0
    ax.tick_params(axis='y', which='major', pad=0)

    #ax.yaxis.set_minor_formatter(ticker.ScalarFormatter()) 
    #ax.yaxis.set_minor_formatter(ticker.NullFormatter()) 
    #ax.yaxis.set_major_formatter(ticker.ScalarFormatter()) 

    #ax.xaxis.set_major_formatter(ticker.ScalarFormatter()) 
    #ax.xaxis.set_minor_formatter(ticker.ScalarFormatter()) 
    '''

    fig.set_size_inches(3,1.8)
    fig.savefig(hp.plot_dir() + name + ".png", bbox_inches = 'tight', dpi = 300)
    fig.savefig(hp.plot_dir() + name + ".pdf", bbox_inches = 'tight', pad_inches = 0, dpi = 300)
    plt.close()

ndiscs = [10,20,30]
version = 1
names = ["SCBA_ndisc" + str(ndisc) + "_v" + str(version) for ndisc in ndiscs]

for name in names:
    plot(name)
