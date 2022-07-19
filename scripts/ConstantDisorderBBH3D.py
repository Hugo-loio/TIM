import helper as hp
import matplotlib.pyplot as plt
import numpy as np

def plot(name):
    fig, ax = plt.subplots()
    data = hp.readfile(name + ".dat")

    data = data[:, data[0,:].argsort()]
    qyz = data[1::3]
    qxz = data[2::3]
    qxy = data[3::3]
    q = 8*np.absolute(np.multiply(np.multiply(qyz,qxz),qxy))

    ax.errorbar(data[0], np.average(q,axis=0), yerr = np.std(q, axis = 0)/np.sqrt(len(q)), capsize = 5, linestyle='-')

    ax.set(xlabel = r'$L$', ylabel = r'$Q$')

    fig.savefig(hp.plot_dir() + name + ".png", dpi = 200)
    plt.show()
    plt.close()

plot("constantDisorderBBH3Dquad_intra1.1_w3")
