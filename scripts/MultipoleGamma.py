import helper as hp
import matplotlib.pyplot as plt

def MultipoleGamma(fname, show: bool, ylab):
    data = hp.readfile(fname)

    fig, ax = plt.subplots()
    for i in range(1,len(data)):
        ax.plot(data[0],data[i],'.')

    ax.set(xlabel = r'$\frac{\gamma}{\lambda}$', ylabel = ylab)

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 200)
    if(show):
        plt.show()
    plt.close()

#quadrupole = ["QuadrupoleBBH2D_20_20","QuadrupoleBBH2D_10_10", "QuadrupoleBBH2DSupercell_4_4", "QuadrupoleBBH2DSupercellWrong_4_4"] 
quadrupole = ["QuadrupoleManyBodyBBH2D_10x10", "QuadrupoleManyBodyBBH2D_5x5", "QuadrupoleManyBodyBBH2D_15x15", "QuadrupoleManyBodyBBH2D_20x20"] 

octupole = ["OctupoleManyBodyBBH3D_4x4x4", "OctupoleManyBodyBBH3D_5x5x5", "OctupoleManyBodyBBH3D_6x6x6"]

boundaryPolx = ["BoundaryPxBBH2D_50x100"]

boundaryPoly = ["BoundaryPyBBH2D_50x100"]

boundaryQxy = ["BoundaryQxyBBH3D_10x10x10", "BoundaryQxyBBH3D_15x15x15"]

boundaryQxz = ["BoundaryQxzBBH3D_10x10x10", "BoundaryQxzBBH3D_15x15x15"]

boundaryQyz = ["BoundaryQyzBBH3D_10x10x10", "BoundaryQyzBBH3D_15x15x15"]

for name in quadrupole:
    MultipoleGamma(name + ".dat", False, r'$q_{xy}$')

for name in octupole:
    MultipoleGamma(name + ".dat", False, r'$o_{xyz}$')

for name in boundaryPolx:
    MultipoleGamma(name + ".dat", False, r'$p_x$')

for name in boundaryPoly:
    MultipoleGamma(name + ".dat", False, r'$p_y$')

for name in boundaryQxy:
    MultipoleGamma(name + ".dat", False, r'$q_{xy}$')

for name in boundaryQxz:
    MultipoleGamma(name + ".dat", False, r'$q_{xz}$')

for name in boundaryQyz:
    MultipoleGamma(name + ".dat", False, r'$q_{yz}$')
