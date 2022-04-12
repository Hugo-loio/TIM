import helper as hp
import matplotlib.pyplot as plt

def WannierBands2D(fname, show: bool, xlab: str):
    data = hp.readfile(fname)

    fig, ax = plt.subplots()
    for i in range(1,len(data)):
        ax.plot(data[0],data[i])

    ax.set(xlabel = xlab, ylabel = r'$\nu$')

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 200)
    if(show):
        plt.show()

names_x = ["WannierBandsBBH2D_inter3_x", "WannierBandsBBH2D_inter2_x", "WannierBandsBBH2D_inter1.5_x", "WannierBandsBBH2D_inter1_x", "WannierBandsBBH2D_inter0.5_x"] 
names_y = ["WannierBandsBBH2D_inter3_y", "WannierBandsBBH2D_inter2_y", "WannierBandsBBH2D_inter1.5_y", "WannierBandsBBH2D_inter1_y", "WannierBandsBBH2D_inter0.5_y"] 

for name in names_x:
    WannierBands2D(name + ".dat", True, r'$k_y$')

for name in names_y:
    WannierBands2D(name + ".dat", True, r'$k_x$')

