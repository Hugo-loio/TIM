import helper as hp
import matplotlib.pyplot as plt

def EnergyBands2D(fname, show: bool):
    data = hp.readfile(fname)

    fig, ax = plt.subplots()
    for i in range(1,len(data)):
        ax.plot(data[0],data[i])

    ax.set(xlabel = r'k', ylabel = r'E', title=r'Energy bands')

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 200)
    if(show):
        plt.show()

names = ["EnergyBandsSSH_inter2", "EnergyBandsSSH_inter1", "EnergyBandsSSH_inter0.5"] 

for name in names:
    EnergyBands2D(name + ".dat", True)

