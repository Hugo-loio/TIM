import helper as hp
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

def matrixHist(fname, show: bool, size: int, bsize: int):
    data = hp.readfile(fname)

    fig, ax = plt.subplots()
    h = ax.hist2d(data[0], data[1], bins=int(size/bsize), range=[[0,size],[0,size]], weights=data[2]) 
    ax.set(xlabel = r'col', ylabel = r'line')
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')

    col = fig.colorbar(h[3], label = 'Weight')

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 200, bbox_inches='tight')
    #ax.colorbar()
    if(show):
        plt.show()

names = ["testH"]

for name in names:
    matrixHist(name + ".dat", True, 100, 20)
