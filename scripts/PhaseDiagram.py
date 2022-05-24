import helper as hp
import matplotlib.pyplot as plt
import numpy as np

def average(data):
    np.average(data[1:], axis=0)

plot_name = ""
names = []
labels = []


fig, ax = plt.subplots()

'''
for name in names:
    data = hp.readfile(fname)

    for i in range(1,len(data)):
        ax.plot(data[0],data[i],'.')

ax.set(xlabel = r'$\frac{\gamma}{\lambda}$', ylabel = ylab)

fig.savefig(hp.plot_dir() + name + ".png", dpi = 200)
if(show):
    plt.show()
    '''
