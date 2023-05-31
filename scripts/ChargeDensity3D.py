import helper as hp
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

plt.rcParams.update({
    "text.usetex": True,
    #"font.family": "Helvetica"
    #"font.size" : 16
})

def ChargeDensity3D(fname, show: bool):
    data = hp.readfile(fname)

    fig = plt.figure()
    ax = fig.add_subplot(projection = '3d')
    img = ax.scatter(data[0],data[1],data[2], c = data[3] , cmap = cm.coolwarm)

    ax.set(xlabel = r'$x$', ylabel = r'$y$', zlabel = r'$z$')

    axis_vecx = np.arange(1, data[0][-1] + 1, 1)
    #axis_vecx[0] += 1
    ax.set_xticks(axis_vecx)
    ax.set_yticks(axis_vecx)
    ax.set_zticks(axis_vecx)
    #ax.set_xticks(np.arange(0, data[0][-1] + 1, 2))
    #ax.set_yticks(np.arange(0, data[1][-1] + 1, 2))
    #ax.set_zticks(np.arange(0, data[2][-1] + 1, 2))

    ax.set_xlabel(r'$x$', labelpad = -4)
    ax.set_ylabel(r'$y$', labelpad = -4)
    ax.set_zlabel(r'$z$', labelpad = -6)

    ax.tick_params(axis='both', which='major', pad=-4)
    ax.tick_params(axis='z', which='major', pad=-2)

    cbar = fig.colorbar(img, shrink = 0.5, pad = -0.01)
    cbar.set_label("Electron Density")
    cbar.ax.ticklabel_format(useOffset=False)
    #pos = cbar.get_position()

    ax.view_init(10,30)
    fig.savefig(hp.plot_dir() + name + ".png", dpi = 200, bbox_inches='tight')
    fig.savefig(hp.plot_dir() + name + ".pdf", bbox_inches = 'tight', pad_inches = 0.2)
    if(show):
        plt.show()

#names = ["ChargeDensityBBH3D_inter2", "ChargeDensityBBH3D_inter1", "ChargeDensityBBH3D_inter0.5","ChargeDensityBBH3D_inter2_delta", "ChargeDensityBBH3D_inter1_delta", "ChargeDensityBBH3D_inter0.5_delta"] 
names = ["ChargeDensityBBH3D_inter2_delta"] 

for name in names:
    ChargeDensity3D(name + ".dat", False)
