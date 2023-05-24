import helper as hp
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": True
    #"font.family": "Helvetica"
    })

def EnergyBands3D(fname, show: bool):
    data = hp.readfile(fname)

    fig = plt.figure()
    ax = fig.add_subplot(projection = '3d')
    for i in range(2,len(data)):
        ax.plot_trisurf(data[0],data[1],data[i], rasterized = True)

    ax.set(xlabel = r'$k_x$', ylabel = r'$k_y$', zlabel = r'$E$')
    ax.set_zlabel(r'$E$', labelpad = -7)
    ax.tick_params(axis='both', which='major', pad=-3)

    ax.view_init(10,-30)
    fig.savefig(hp.plot_dir() + name + ".png", dpi = 300, bbox_inches='tight', pad_inches = 0.2)
    fig.savefig(hp.plot_dir() + name + ".pdf", dpi = 300, bbox_inches = 'tight', pad_inches = 0.2)
    if(show):
        plt.show()

#names = ["EnergyBandsSOTAI_m2", "EnergyBandsSOTAI_m1", "EnergyBandsSOTAI_m0.5", "EnergyBandsSOTAI_m1.01", "EnergyBandsSOTAI_m0.99"]
#names = ["EnergyBandsBBH2D_inter2", "EnergyBandsBBH2D_inter1.5", "EnergyBandsBBH2D_inter1", "EnergyBandsBBH2D_inter0.5", "EnergyBandsSSH2D_inter2", "EnergyBandsSSH2D_inter1", "EnergyBandsSSH2D_inter0.5"] 
names = ["EnergyBandsBBH2D_inter2", "EnergyBandsBBH2D_inter1.5", "EnergyBandsBBH2D_inter1", "EnergyBandsBBH2D_inter0.5"]

for name in names:
    EnergyBands3D(name + ".dat", False)
