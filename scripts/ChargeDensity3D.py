import helper as hp
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

def ChargeDensity3D(fname, show: bool):
    data = hp.readfile(fname)

    fig = plt.figure()
    ax = fig.add_subplot(projection = '3d')
    img = ax.scatter(data[0],data[1],data[2], c = data[3] , cmap = cm.coolwarm)

    ax.set(xlabel = r'$x$', ylabel = r'$y$', zlabel = r'$z$')

    ax.set_xticks(np.arange(0, data[0][-1] + 1, 2))
    ax.set_yticks(np.arange(0, data[1][-1] + 1, 2))
    ax.set_zticks(np.arange(0, data[2][-1] + 1, 2))

    cbar = fig.colorbar(img, shrink = 0.4)
    cbar.set_label("Electron Density")
    cbar.ax.ticklabel_format(useOffset=False)

    ax.view_init(10,30)
    fig.savefig(hp.plot_dir() + name + ".png", dpi = 200, bbox_inches='tight')
    if(show):
        plt.show()

names = ["ChargeDensityBBH3D_inter2", "ChargeDensityBBH3D_inter1", "ChargeDensityBBH3D_inter0.5","ChargeDensityBBH3D_inter2_delta", "ChargeDensityBBH3D_inter1_delta", "ChargeDensityBBH3D_inter0.5_delta"] 

for name in names:
    ChargeDensity3D(name + ".dat", False)
