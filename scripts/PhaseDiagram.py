import helper as hp
import matplotlib.pyplot as plt
import numpy as np

def plotPol():
    data = hp.readfile(namesPol[i] + ".dat")

    data = data[:, data[0,:].argsort()]
    px = data[1:-2:2]
    py = data[2:-1:2]
    p = np.absolute(np.multiply(px,py))

    ax.errorbar(data[0], np.average(px,axis=0), yerr = np.std(px,axis=0)/np.sqrt(np.size(data)), label = r'$p_x$, ' + labelsQuad[i], capsize = 5, linestyle='-')
    ax.errorbar(data[0], np.average(py,axis=0), yerr = np.std(py,axis=0)/np.sqrt(np.size(data)), label = r'$p_y$, ' + labelsQuad[i], capsize = 5, linestyle='-')
    ax.errorbar(data[0], np.average(p,axis=0), yerr = np.std(p,axis=0)/np.sqrt(np.size(data)), label = r'$P$, ' + labelsQuad[i], capsize = 5, linestyle='-')

def plotQuad():
    data = hp.readfile(namesQuad[i] + ".dat")

    data = data[:, data[0,:].argsort()]
    ax.errorbar(data[0], np.average(data[1:],axis=0), yerr = np.std(data[1:],axis=0)/np.sqrt(np.size(data)), label = labelsQuad[i], capsize = 5, linestyle='-')


plot_name = "phaseDiagramSOTAI"
#namesQuad = ["phaseDiagramSOTAI_10x10_m1.1", "phaseDiagramSOTAI_20x20_m1.1"]
namesQuad = ["phaseDiagramSOTAI_10x10_m1.1"]
labelsQuad = ["L = 10", "L = 20"]
namesPol = ["phaseDiagramSOTAIpol_10x10_m1.1"]
labelsPol = ["L = 10"]


fig, ax = plt.subplots()

'''
for i in range(0, len(namesQuad)):
    plotQuad()
'''

for i in range(0, len(namesPol)):
    plotPol()

ax.set(xlabel = r'$W$')
ax.legend()

fig.savefig(hp.plot_dir() + plot_name + ".png", dpi = 200)
plt.show()
