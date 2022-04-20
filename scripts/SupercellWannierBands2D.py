import helper as hp
import matplotlib.pyplot as plt

def SupercellWannierBands(fname, show: bool, xlab: str):
    data = hp.readfile(fname)

    fig, ax = plt.subplots()
    for i in range(1,len(data)):
        if(i < ((len(data)-1)/2 + 1)):
            ax.plot(data[0],data[i],color = 'r')
        else:
            ax.plot(data[0],data[i],color = 'b')

    ax.set(xlabel = xlab, ylabel = r'$\nu$')

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 200)
    if(show):
        plt.show()

names_x = ["SupercellWannierBandsBBH2D_inter2_x", "SupercellWannierBandsBBH2D_inter1_x", "SupercellWannierBandsBBH2D_inter0.5_x"] 
names_y = ["SupercellWannierBandsBBH2D_inter2_y"] 


for name in names_x:
    SupercellWannierBands(name + ".dat", True, r'$\theta_y$')

for name in names_y:
    SupercellWannierBands(name + ".dat", False, r'$\theta_x$')

