import helper as hp
import matplotlib.pyplot as plt
import numpy as np
import ConstantDisorderBBH3D as extra

plt.style.use('science')

def plotQuad(name, label, ax, detail):
    data = hp.readfile(name + ".dat")

    data = data[:, data[0,:].argsort()]
    qyz = data[1::3]
    qxz = data[2::3]
    qxy = data[3::3]
    q = 8*np.absolute(np.multiply(np.multiply(qyz,qxz),qxy))

    if(detail):
        ax.plot(data[0], np.average(qyz,axis=0), label = r'$q_{yz}$, ' + label,  linestyle='-')
        ax.plot(data[0], np.average(qxz,axis=0), label = r'$q_{xz}$, ' + label,  linestyle='-')
        ax.plot(data[0], np.average(qxy,axis=0), label = r'$q_{xy}$, ' + label,  linestyle='-')

    ax.plot(data[0], np.average(q,axis=0), label = label,  linestyle='-')

def plot(name, fileNames, labels, detail = True, show = True):
    fig, ax = plt.subplots()

    for i in range(0, len(fileNames)):
        plotQuad(fileNames[i], labels[i], ax, detail)

    ax.errorbar(extra.res[0], extra.res[1], yerr = extra.res[2], label = 'extrapolation', linestyle='', marker = ".", markersize=5, markeredgewidth = 1, linewidth = 1, capsize = 2, capthick = 1, color = 'darkorchid')
    #print(extra.res)

    ax.set(xlabel = r'$W$', ylabel = r'$Q$')
    ymin, ymax = ax.get_ylim()
    #ax.margins(x = 0)
    #ax.set_ylim([ymin, ymax*1.2])
    hp.totaiPhases(ax)
    ax.legend(bbox_to_anchor=(0.5,0.4), fontsize=8)

    fig.set_size_inches(3.4, 2.5)
    fig.savefig(hp.plot_dir() + name + ".png", dpi = 300, bbox_inches = 'tight', pad_inches = 0.01)
    fig.savefig(hp.plot_dir() + name + ".pdf", bbox_inches = 'tight', pad_inches = 0.01)
    if(show):
        plt.show()
    #plt.close()


namesQuad = ["L5", "L7", "L10", "L12", "L13", "L15", "L17", "L18", "L20", "L22", "L24"]
namesQuad = ["phaseDiagramBBH3Dquad_" + name + "_intra1.1" for name in namesQuad]
labelsQuad = ["L = 5", "L = 7", "L = 10", "L = 12", "L = 13", "L = 15", "L = 17", "L = 18", "L = 20", "L = 22", "L = 24"]

#plot("PhaseDiagramBBH3D_intra1.1", namesQuad[2:], labelsQuad[2:], True, False)
#plot("PhaseDiagramBBH3D_intra1.1_Q", namesQuad[2:], labelsQuad[2:], False, False)

sizesQuad = ["18", "20", "22", "24"]
namesQuad = ["phaseDiagramBBH3Dquad_L" + name + "_intra1.1" for name in sizesQuad]
labelsQuad = ["L = " + name for name in sizesQuad]
plot("PhaseDiagramBBH3Dpretty", namesQuad, labelsQuad, False, False)

namesQuad2 = ["phaseDiagramBBH3Dquad_L7_intra0.9"]
labelsQuad2 = ["L=7"]

#plot("PhaseDiagramBBH3D_intra0.9", namesQuad2, labelsQuad2, True, False)

namesQuad3 = ["phaseDiagramBBH3Dquad_L7_intra0.5"]
labelsQuad3 = ["L=7"]

#plot("PhaseDiagramBBH3D_intra0.5", namesQuad3, labelsQuad3, True, False)

namesQuad4 =["", "_twistPi"]
namesQuad4 =["phaseDiagramBBH3Dquad_L15_intra1.1" + name for name in namesQuad4]
labelsQuad4 =[r'$L = 15, \theta = 0$', r'$L = 15, \theta = \pi$']
#plot("PhaseDiagramBBH3D_L15_intra1.1_twist", namesQuad4, labelsQuad4, False, False)
