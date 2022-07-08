import helper as hp
import matplotlib.pyplot as plt
import numpy as np

def plotQuad(name, label, ax, detail):
    data = hp.readfile(name + ".dat")

    data = data[:, data[0,:].argsort()]
    qyz = data[1:-3:3]
    qxz = data[2:-2:3]
    qxy = data[3:-1:3]
    q = 8*np.absolute(np.multiply(np.multiply(qyz,qxz),qxy))

    if(detail):
        ax.plot(data[0], np.average(qyz,axis=0), label = r'$q_{yz}$, ' + label,  linestyle='-')
        ax.plot(data[0], np.average(qxz,axis=0), label = r'$q_{xz}$, ' + label,  linestyle='-')
        ax.plot(data[0], np.average(qxy,axis=0), label = r'$q_{xy}$, ' + label,  linestyle='-')
    
    ax.plot(data[0], np.average(q,axis=0), label = r'$Q$, ' + label,  linestyle='-')

def plot(name, fileNames, labels, detail = True):
    fig, ax = plt.subplots()

    for i in range(0, len(fileNames)):
        plotQuad(fileNames[i], labels[i], ax, detail)

    ax.set(xlabel = r'$W$')
    ax.legend()

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 200)
    plt.show()
    plt.close()


namesQuad = ["phaseDiagramBBH3Dquad_L5_intra1.1","phaseDiagramBBH3Dquad_L7_intra1.1","phaseDiagramBBH3Dquad_L10_intra1.1"]
labelsQuad = ["L = 5","L=7","L=10"]

plot("phaseDiagramBBH3D_intra1.1", namesQuad, labelsQuad)
plot("phaseDiagramBBH3D_intra1.1_Q", namesQuad, labelsQuad, False)

namesQuad2 = ["phaseDiagramBBH3Dquad_L7_intra0.9"]
labelsQuad2 = ["L=7"]

plot("phaseDiagramBBH3D_intra0.9", namesQuad2, labelsQuad2)

namesQuad3 = ["phaseDiagramBBH3Dquad_L7_intra0.5"]
labelsQuad3 = ["L=7"]

plot("phaseDiagramBBH3D_intra0.5", namesQuad3, labelsQuad3)
