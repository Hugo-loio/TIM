import helper as hp
import matplotlib.pyplot as plt

def EnergyBands2D(fname, show: bool):
    data = hp.readfile(fname)

    fig, ax = plt.subplots()
    for i in range(1,len(data)):
        #ax.plot(data[0],data[i],'o', markersize=1)
        ax.plot(data[0],data[i])
        ax.axvline(x=1, linestyle = 'dashed', color='black', linewidth='0.8')
        ax.axvline(x=2, linestyle = 'dashed', color='black', linewidth='0.8')
        ax.axvline(x=3, linestyle = 'dashed', color='black', linewidth='0.8')
        ax.axvline(x=4, linestyle = 'dashed', color='black', linewidth='0.8')
        ax.axvline(x=5, linestyle = 'dashed', color='black', linewidth='0.8')
        ax.axvline(x=6, linestyle = 'dashed', color='black', linewidth='0.8')
        plt.xlim([0,6])

    ax.set(xlabel = '', ylabel = r'E', title=r'')

    #labels = [item.get_text() for item in ax.get_xticklabels()]
    labels = [r'$\Gamma$', r'$X$', r'$M$', r'$\Gamma$', r'$R$', r'$X|M$', r'$R$']
    ax.set_xticklabels(labels)

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 200, bbox_inches='tight')
    if(show):
        plt.show()

names = ["EnergyBandsBBH3D_inter2","EnergyBandsBBH3D_inter1","EnergyBandsBBH3D_inter0.5"] 

for name in names:
    EnergyBands2D(name + ".dat", True)
