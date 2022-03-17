import helper as hp
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

def ChargeDensity2D(fname, show: bool):
    data = hp.readfile(fname)

    fig = plt.figure()
    ax = fig.add_subplot(projection = '3d')
    for i in range(2,len(data)):
        ax.plot_trisurf(data[0],data[1],data[i], cmap = cm.coolwarm)

    ax.set(xlabel = r'$x$', ylabel = r'$y$', zlabel = r'Electron density')

    ax.set_xticks(np.arange(0, data[0][-1] + 1, 5))
    ax.set_yticks(np.arange(0, data[1][-1] + 1, 5))

    ax.view_init(10,30)
    fig.savefig(hp.plot_dir() + name + ".png", dpi = 200, bbox_inches='tight')
    if(show):
        plt.show()

names = ["ChargeDensityBBH2D_inter2", "ChargeDensityBBH2D_inter1", "ChargeDensityBBH2D_inter0.5","ChargeDensityBBH2D_inter2_delta", "ChargeDensityBBH2D_inter1_delta", "ChargeDensityBBH2D_inter0.5_delta"] 

for name in names:
    ChargeDensity2D(name + ".dat", True)
