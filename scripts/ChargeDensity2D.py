import helper as hp
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

plt.rcParams.update({
    "text.usetex": True,
    #"font.family": "Helvetica"
    #"font.size" : 16
})

def ChargeDensity2D(fname, show: bool):
    data = hp.readfile(fname)

    fig = plt.figure()
    ax = fig.add_subplot(projection = '3d')
    for i in range(2,len(data)):
        ax.plot_trisurf(data[0],data[1],data[i], cmap = cm.coolwarm, vmin = 1.85, vmax = 2.15)

    ax.view_init(10,30)
    ax.set(xlabel = r'$x$', ylabel = r'$y$', zlabel = r'Electron density')

    axis_vecx = np.arange(0, data[0][-1] + 1, 5)
    axis_vecx[0] += 1
    ax.set_xticks(axis_vecx)
    axis_vecy = np.arange(0, data[1][-1] + 1, 5)
    axis_vecy[0] += 1
    ax.set_yticks(axis_vecy)
    ax.ticklabel_format(useOffset=False)
    ax.set_zlim(1.7,2.3)

    ax.tick_params(axis='both', which='major', pad=-4)
    ax.tick_params(axis='z', which='major', pad=-2)
    ax.set_xlabel(r'$x$', labelpad = -4)
    ax.set_ylabel(r'$y$', labelpad = -4)
    ax.set_zlabel(r'Electron density', labelpad = -6)

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 300, bbox_inches='tight', pad_inches = 0.3)
    fig.savefig(hp.plot_dir() + name + ".pdf", bbox_inches = 'tight', pad_inches = 0.2)
    if(show):
        plt.show()

#names = ["ChargeDensityBBH2D_inter2", "ChargeDensityBBH2D_inter1", "ChargeDensityBBH2D_inter0.5","ChargeDensityBBH2D_inter2_delta", "ChargeDensityBBH2D_inter1_delta", "ChargeDensityBBH2D_inter0.5_delta", "ChargeDensityBBH2D_inter1.01_delta", "ChargeDensityBBH2D_inter0.99_delta"] 

names = ["ChargeDensityBBH2D_inter2_delta"] 

#names = ["ChargeDensitySOTAI_m0.5", "ChargeDensitySOTAI_m2", "ChargeDensitySOTAI_m0.99", "ChargeDensitySOTAI_m1.01", "ChargeDensitySOTAI_m1.1_w0", "ChargeDensitySOTAI_m1.1_w1", "ChargeDensitySOTAI_m1.1_w2.5", "ChargeDensitySOTAI_m1.1_w3","ChargeDensitySOTAI_m0.5_w0", "ChargeDensitySOTAI_m0.5_w1", "ChargeDensitySOTAI_m0.5_w2", "ChargeDensitySOTAI_m0.5_w3"]

for name in names:
    ChargeDensity2D(name + ".dat", False)
