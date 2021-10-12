import helper as hp
import matplotlib.pyplot as plt

name = "1DSSHEBands"
data = hp.readfile(name + ".dat")

fig, ax = plt.subplots()
ax.plot(data[0],data[1], color = 'blue')
ax.plot(data[0],data[2], color = 'orange')

ax.set(xlabel = r'k', ylabel = r'E', title=r'SSH Energy bands')

fig.savefig(hp.plot_dir() + name + ".png", dpi = 200)
plt.show()
