import helper as hp
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

plt.style.use('science')

def ldos(fname, show: bool):
    data = hp.readfile(fname)
    x = data[0::3,0]
    y = data[1::3,0]
    zVec = np.average(data[2::3], axis = 1)
    #print(zVec)
    #print(z)
    zMat = np.empty([int(x[-1])+1, int(y[-1])+1])
    zMin = 0
    zMax = 1
    for i in range(len(x)):
        #print([i,j])
        z = zVec[int(y[i])*int(x[-1]) + int(x[i])]
        zMat[int(x[i]),int(y[i])] = z
        if(z < zMin):
            zMin = z
        elif(z > zMax):
            zMax = z

    fig, ax = plt.subplots();
    im = ax.imshow(zMat, vmin = zMin, vmax = zMax)
    fig.colorbar(im , ax = ax)

    ax.set(xlabel = r'$x$', ylabel = r'$y$')

    ax.set_xticks(np.arange(0, data[0][-1] + 1, 5))
    ax.set_yticks(np.arange(0, data[1][-1] + 1, 5))
    ax.ticklabel_format(useOffset=False)
    #ax.set_zlim(1.7,2.3)

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 300, bbox_inches='tight')
    if(show):
        plt.show()

weights = ['1', '2.6', '3.4']
names = ['ldosSOTAI_L10_w' + w + '_E0_nMu1024_m1.1' for w in weights]

for name in names:
    ldos(name + ".dat", False)
