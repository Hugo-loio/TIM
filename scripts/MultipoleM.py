import helper as hp
import matplotlib.pyplot as plt

def MultipoleM(fname, show: bool, ylab):
    data = hp.readfile(fname)

    fig, ax = plt.subplots()
    for i in range(1,len(data)):
        ax.plot(data[0],data[i],'.')

    ax.set(xlabel = r'$m$', ylabel = ylab)

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 200)
    if(show):
        plt.show()

#quadrupole = ["QuadrupoleSOTAI_10_10", "QuadrupoleSOTAISupercell_4_4"] 
#quadrupole = ["QuadrupoleManyBodySOTAI_10x10","QuadrupoleManyBodySOTAI_20x20"] 
quadrupole = ["QuadrupoleManyBodyDisSOTAI_10x10","QuadrupoleManyBodyDisSOTAI_20x20"] 

octupole = []

for name in quadrupole:
    MultipoleM(name + ".dat", False, r'$q_{xy}$')

for name in octupole:
    MultipoleM(name + ".dat", False, r'$o_{xyz}$')

