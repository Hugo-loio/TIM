import helper as hp
import matplotlib.pyplot as plt
import numpy as np

plt.style.use('science')

def plotPol():
    data = hp.readfile(namesPol[i] + ".dat")

    data = data[:, data[0,:].argsort()]
    px = data[1::2]
    py = data[2::2]
    p = 4*np.absolute(np.multiply(px,py))

    #ax.errorbar(data[0], np.average(px,axis=0), yerr = np.std(px,axis=0)/np.sqrt(np.size(data)), label = r'$p_x$, ' + labelsPol[i], capsize = 5, linestyle='-')
    #ax.errorbar(data[0], np.average(px,axis=0), yerr = np.std(px,axis=0)/np.sqrt(np.size(data)), label = r'$p_x$, ' + labelsPol[i], capsize = 5, linestyle='-')
    #ax.errorbar(data[0], np.average(py,axis=0), yerr = np.std(py,axis=0)/np.sqrt(np.size(data)), label = r'$p_y$, ' + labelsPol[i], capsize = 5, linestyle='-')
    ax.errorbar(data[0], np.average(p,axis=0), label = r'$P$, ' + labelsPol[i], capsize = 5, linestyle='-')

def plotQuad():
    data = hp.readfile(namesQuad[i] + ".dat")

    data = data[:, data[0,:].argsort()]
    ax.errorbar(data[0], np.average(data[1:],axis=0), label = r'$q_{xy}$, ' + labelsQuad[i], capsize = 5, linestyle='-')


plot_name = "PhaseDiagramSOTAI"
#namesQuad = ["phaseDiagramSOTAI_10x10_m1.1", "phaseDiagramSOTAI_20x20_m1.1"]
namesQuad = ["phaseDiagramSOTAI_20x20_m1.1","phaseDiagramSOTAI_40x40_m1.1","phaseDiagramSOTAI_60x60_m1.1"]
labelsQuad = ["L = 20", "L = 40", "L = 60"]
#namesPol = ["phaseDiagramSOTAIpol_10x10_m1.1", "phaseDiagramSOTAIpol_50x50_m1.1"]
#labelsPol = ["L = 10", "L = 50"]
namesPol = ["phaseDiagramSOTAIpol_300x300_m1.1", "phaseDiagramSOTAIpol_400x400_m1.1", "phaseDiagramSOTAIpol_500x500_m1.1"]
labelsPol = ["L = 300", "L = 400", "L = 500"]


fig, ax = plt.subplots()

for i in range(0, len(namesQuad)):
    plotQuad()

for i in range(0, len(namesPol)):
    plotPol()

ax.set(xlabel = r'$W$')
ax.margins(x = 0)
ax.set_ylim(-0.04, 1.15)
hp.sotaiPhases(ax)
ax.legend(fontsize = 8, bbox_to_anchor=(0.66,0.8))

fig.set_size_inches(4, 3)
fig.savefig(hp.plot_dir() + plot_name + ".png", bbox_inches = 'tight', dpi = 300)
fig.savefig(hp.plot_dir() + plot_name + ".pdf", bbox_inches = 'tight', pad_inches = 0.01)
#plt.show()
