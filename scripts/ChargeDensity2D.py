import helper as hp
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

def ChargeDensity2D(fname, show: bool):
    data = hp.readfile(fname)

    fig = plt.figure()
    ax = fig.add_subplot(projection = '3d')
    for i in range(2,len(data)):
        ax.plot_trisurf(data[0],data[1],data[i], cmap = cm.coolwarm, vmin = 1.85, vmax = 2.15)

    ax.set(xlabel = r'$x$', ylabel = r'$y$', zlabel = r'Electron density')

    ax.set_xticks(np.arange(0, data[0][-1] + 1, 5))
    ax.set_yticks(np.arange(0, data[1][-1] + 1, 5))
    ax.ticklabel_format(useOffset=False)
    ax.set_zlim(1.7,2.3)

    ax.view_init(10,30)
    fig.savefig(hp.plot_dir() + name + ".png", dpi = 200, bbox_inches='tight')
    if(show):
        plt.show()

#names = ["ChargeDensityBBH2D_inter2", "ChargeDensityBBH2D_inter1", "ChargeDensityBBH2D_inter0.5","ChargeDensityBBH2D_inter2_delta", "ChargeDensityBBH2D_inter1_delta", "ChargeDensityBBH2D_inter0.5_delta", "ChargeDensityBBH2D_inter1.01_delta", "ChargeDensityBBH2D_inter0.99_delta"] 

names = ["ChargeDensitySOTAI_m0.5", "ChargeDensitySOTAI_m2", "ChargeDensitySOTAI_m0.99", "ChargeDensitySOTAI_m1.01", "ChargeDensitySOTAI_m1.1_w0", "ChargeDensitySOTAI_m1.1_w1", "ChargeDensitySOTAI_m1.1_w2.5", "ChargeDensitySOTAI_m1.1_w3","ChargeDensitySOTAI_m0.5_w0", "ChargeDensitySOTAI_m0.5_w1", "ChargeDensitySOTAI_m0.5_w2", "ChargeDensitySOTAI_m0.5_w3"]

for name in names:
    ChargeDensity2D(name + ".dat", False)
