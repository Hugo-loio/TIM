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
    plt.close()

names = ["testH"]

names2 = ["testH2"]

names1 = ["hamBBH3D_4x4x4_noPBC", "hamBBH3D_4x4x4_PBCxy", "hamBBH3D_4x4x4_PBCx", "hamBBH3D_4x4x4_PBCy"]

boundaryHamBBH2D = ["BoundaryHamxBBH2D_10x10_intra0.5", "BoundaryHamyBBH2D_10x10_intra0.5","BoundaryHamxBBH2D_10x10_intra2", "BoundaryHamyBBH2D_10x10_intra2", "BoundaryHamxBBH2D_10x5_intra0.5","BoundaryHamxBBH2D_10x1_intra0.5","BoundaryHamxBBH2D_10x3_intra0.5", "BoundaryHamxBBH2D_10x5_intra2","BoundaryHamxBBH2D_10x1_intra2","BoundaryHamxBBH2D_10x3_intra2"]

'''
for name in names:
    matrixHist(name + ".dat", True, 512, 64)

for name in names2:
    matrixHist(name + ".dat", False, 16, 1)

for name in names1:
    matrixHist(name + ".dat", False, 512, 8)
'''

for name in boundaryHamBBH2D:
    matrixHist(name + ".dat", False, 20, 1)
