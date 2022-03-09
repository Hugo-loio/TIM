import helper as hp
import matplotlib.pyplot as plt

def WannierBands2D(fname, show: bool):
    data = hp.readfile(fname)

    fig, ax = plt.subplots()
    for i in range(1,len(data)):
        ax.plot(data[0],data[i])

    ax.set(xlabel = r'k', ylabel = r'E', title=r'Wannier Bands')

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 200)
    if(show):
        plt.show()

names = ["WannierBandsBBH2D_inter2_x", "WannierBandsBBH2D_inter2_y","WannierBandsBBH2D_inter0.5_x", "WannierBandsBBH2D_inter0.5_y"] 

for name in names:
    WannierBands2D(name + ".dat", True)

