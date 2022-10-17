import helper as hp
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

plt.style.use('science')

def plot(name, fileName, show = True):
    data = hp.readfile(fileName + ".dat")
    nPointsW = len(np.where(data[1] == data[1][0])[0])
    nPointsM = len(np.where(data[0] == data[0][0])[0])
    print(nPointsM)
    print(nPointsW)
    fig, ax = plt.subplots()

    colormap = plt.cm.nipy_spectral
    colors = [colormap(i) for i in np.linspace(0, 1, nPointsW)]
    ax.set_prop_cycle('color', colors)

    for m in range(nPointsW):
        dataPlot = data[1:, np.where(data[0] == data[0][m*nPointsM])[0]]
        ax.plot(dataPlot[0], dataPlot[1], label = str(data[0][m*nPointsM]), linestyle = '-')

    ax.set(xlabel = "L", ylabel = r'$\Lambda$')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([1, 14])
    ax.legend(fontsize = 4, ncol = 1, loc = 'upper right', title = r'$W$', title_fontsize = 6)
    #ax.xaxis.set_minor_formatter(mticker.ScalarFormatter())
    #ax.yaxis.set_minor_formatter(mticker.ScalarFormatter())

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 300)
    fig.savefig(hp.plot_dir() + name + ".eps")
    if(show):
        plt.show()
    plt.close()


plot("Anderson3DTMM", "tmmAnderson3D_E0_t1", False)
