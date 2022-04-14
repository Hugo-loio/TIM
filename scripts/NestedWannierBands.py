import helper as hp
import matplotlib.pyplot as plt

def NestedWannierBands(fname, show: bool, xlab: str):
    data = hp.readfile(fname)

    fig, ax = plt.subplots()
    for i in range(1,len(data)):
        ax.plot(data[0],data[i])

    ax.set(xlabel = xlab, ylabel = r'$\nu$')

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 200)
    if(show):
        plt.show()

names_z = ["NestedWannierBandsBBH3D_inter2_z","NestedWannierBandsBBH3D_inter1_z","NestedWannierBandsBBH3D_inter0.5_z"] 
names_y = [] 


for name in names_z:
    NestedWannierBands(name + ".dat", True, r'$k_z$')

for name in names_y:
    NestedWannierBands(name + ".dat", False, r'$k_y$')

