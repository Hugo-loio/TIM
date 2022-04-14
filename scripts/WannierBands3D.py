import helper as hp
import matplotlib.pyplot as plt

def WannierBands3D(fname, show: bool, xlab: str, ylab: str):
    data = hp.readfile(fname)

    fig = plt.figure()
    ax = fig.add_subplot(projection = '3d')
    for i in range(2,len(data)):
        ax.plot_trisurf(data[0],data[1],data[i])

    ax.set(xlabel = xlab, ylabel = ylab, zlabel = r'$\nu$')

    ax.view_init(10,30)
    fig.savefig(hp.plot_dir() + name + ".png", dpi = 200, bbox_inches='tight')
    if(show):
        plt.show()

names_x = ["WannierBandsBBH3D_inter2_x", "WannierBandsBBH3D_inter1.01_x", "WannierBandsBBH3D_inter1_x", "WannierBandsBBH3D_inter0.99_x", "WannierBandsBBH3D_inter0.5_x"] 
names_y = ["WannierBandsBBH3D_inter2_y"] 
names_z = ["WannierBandsBBH3D_inter2_z"] 

for name in names_x:
    WannierBands3D(name + ".dat", False, r'$k_y$', r'$k_z$')

for name in names_y:
    WannierBands3D(name + ".dat", False, r'$k_x$', r'$k_z$')

for name in names_z:
    WannierBands3D(name + ".dat", False, r'$k_x$', r'$k_y$')
