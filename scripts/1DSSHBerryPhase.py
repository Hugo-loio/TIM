import helper as hp
import matplotlib.pyplot as plt

name = "1DSSHBerryPhase"
data = hp.readfile(name + ".dat")

fig, ax = plt.subplots()
ax.plot(data[0],data[1], color = 'blue')

ax.set(xlabel = r'$J\prime$', ylabel = r'$\phi$', title=r'Berry phase, $J = 1$')

fig.savefig(hp.plot_dir() + name + ".png", dpi = 200)
plt.show()
