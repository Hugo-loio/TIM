import helper as hp
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np


def phasePlot(fname, show: bool, label, part :int = -1):
    data = hp.readfile(fname)
    if(part != -1):
        z = np.average(data[2 + part::3], axis=0)
    else:
        z = np.average(data[2:], axis=0)
    size = int(np.sqrt(len(data[1])))
    x = np.reshape(data[0], (size, size))
    y = np.reshape(data[1], (size, size))
    z = np.reshape(z, (size,size))

    fig, ax = plt.subplots()
    im = ax.pcolormesh(x, y, z, vmin = np.min(z), vmax = np.max(z), cmap=plt.colormaps['cividis'], shading = 'auto')

    ax.set_xlabel(r'W', fontsize = 20)
    ax.set_ylabel(r'$\frac{v}{w}$', fontsize = 20)
    cb = fig.colorbar(im)
    cb.set_label(label, size = 20)

    match part:
        case -1:
            fig.savefig(hp.plot_dir() + name + ".png", dpi = 200, bbox_inches='tight')
        case 0:
            fig.savefig(hp.plot_dir() + name + "_real" + ".png", dpi = 200, bbox_inches='tight')
        case 1:
            fig.savefig(hp.plot_dir() + name + "_orbital" + ".png", dpi = 200, bbox_inches='tight')
        case 2:
            fig.savefig(hp.plot_dir() + name + "_disconnected" + ".png", dpi = 200, bbox_inches='tight')
    if(show):
        plt.show()
    plt.close()

pol = ["L50", "L64", "L100", "L200"]
pol = ["phaseDiagramSSHpol_" + name for name in pol];

ipr = ["L100_n10", "L200_n10", "L500_n10"]
ipr = ["phaseDiagramSSHipr_" + name for name in ipr];

ent = ["L20", "L50", "L64", "L100"]
ent = ["phaseDiagramSSHent_" + name for name in ent];

for name in pol:
    phasePlot(name + ".dat", False, r'P')

for name in ipr:
    phasePlot(name + ".dat", False, r'IPR')

for name in ent:
    phasePlot(name + ".dat", False, r'$S_A/\log2$', 0)
    phasePlot(name + ".dat", False, r'$S_{Orb}$', 1)
    phasePlot(name + ".dat", False, r'$S^D/\log2$', 2)
