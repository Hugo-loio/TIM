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

#quadrupole = ["QuadrupoleBBH2D_20_20","QuadrupoleBBH2D_10_10", "QuadrupoleBBH2DSupercell_4_4", "QuadrupoleBBH2DSupercellWrong_4_4"] 
quadrupole = ["QuadrupoleManyBodyBBH2D_10x10", "QuadrupoleManyBodyBBH2D_5x5", "QuadrupoleManyBodyBBH2D_15x15", "QuadrupoleManyBodyBBH2D_20x20"] 

octupole = ["OctupoleManyBodyBBH3D_4x4x4", "OctupoleManyBodyBBH3D_5x5x5", "OctupoleManyBodyBBH3D_6x6x6"]

for name in quadrupole:
    MultipoleGamma(name + ".dat", False, r'$q_{xy}$')

for name in octupole:
    MultipoleGamma(name + ".dat", False, r'$o_{xyz}$')

