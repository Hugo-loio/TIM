import helper as hp
import matplotlib.pyplot as plt
import numpy as np

def average(data):
    return np.average(data[1:], axis=0)

plot_name = "phaseDiagramSOTAI"
names = ["phaseDiagramSOTAI_10x10_m1.1", "phaseDiagramSOTAI_20x20_m1.1"]
#names = ["phaseDiagramSOTAI_10x10_m1.1"]
labels = ["L = 10", "L = 20"]


fig, ax = plt.subplots()

for i in range(0, len(names)):
    data = hp.readfile(names[i] + ".dat")

    data = data[:, data[0,:].argsort()]
    ax.errorbar(data[0], np.average(data[1:],axis=0), yerr = np.std(data[1:],axis=0)/np.sqrt(np.size(data)), label = labels[i], capsize = 5, linestyle='-')

ax.set(xlabel = r'$W$')
ax.legend()

fig.savefig(hp.plot_dir() + plot_name + ".png", dpi = 200)
plt.show()
