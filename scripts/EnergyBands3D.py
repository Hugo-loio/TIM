import helper as hp
import matplotlib.pyplot as plt

def EnergyBands3D(fname, show: bool):
    data = hp.readfile(fname)

    fig = plt.figure()
    ax = fig.add_subplot(projection = '3d')
    for i in range(2,len(data)):
        ax.plot_trisurf(data[0],data[1],data[i])

    ax.set(xlabel = r'$k_x$', ylabel = r'$k_y$', zlabel = r'$E$', title=r'Energy bands')

    ax.view_init(10,45)
    fig.savefig(hp.plot_dir() + name + ".png", dpi = 200)
    if(show):
        plt.show()

names = ["EnergyBandsBBH2D_inter2", "EnergyBandsBBH2D_inter1", "EnergyBandsBBH2D_inter0.5", "EnergyBandsSSH2D_inter2", "EnergyBandsSSH2D_inter1", "EnergyBandsSSH2D_inter0.5"] 

for name in names:
    EnergyBands3D(name + ".dat", False)
