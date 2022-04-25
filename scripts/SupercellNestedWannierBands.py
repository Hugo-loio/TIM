import helper as hp
import matplotlib.pyplot as plt

def SupercellNestedWannierBands(fname, show: bool, xlab: str):
    data = hp.readfile(fname)

    fig, ax = plt.subplots()
    for i in range(1,len(data)):
        if(i < ((len(data)-1)/2 + 1)):
            ax.plot(data[0],data[i],color = 'r')
        else:
            ax.plot(data[0],data[i],color = 'b')

    ax.set(xlabel = xlab, ylabel = r'$\eta$')

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 200)
    if(show):
        plt.show()

names_z = ["SupercellNestedWannierBandsBBH3D_intra2_1x1x1_z","SupercellNestedWannierBandsBBH3D_intra2_1x1x2_z"] 


for name in names_z:
    SupercellNestedWannierBands(name + ".dat", True, r'$\theta_z$')
