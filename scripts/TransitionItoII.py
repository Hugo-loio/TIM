import helper as hp
import matplotlib.pyplot as plt
import numpy as np
import ConstantDisorderBBH3D as extra

plt.style.use('science')

def getWc(name):
    data = hp.readfile(name + ".dat")

    data = data[:, data[0,:].argsort()]
    qyz = data[1::3]
    qxz = data[2::3]
    qxy = data[3::3]
    q = np.average(8*np.absolute(np.multiply(np.multiply(qyz,qxz),qxy)), axis = 0)
    return data[0][np.argwhere(q > 1E-6)[0][0]]

#sizesQuad = ["12", "20", "22", "24"]
gammas = np.round(np.arange(1.025, 1.285, 0.025),3)
gammas = [str(gamma) for gamma in gammas]

namesQuad = ["phaseDiagramBBH3Dquad_L12_intra" + gamma for gamma in gammas]
for i,name in enumerate(namesQuad):
    Wc = getWc(name)
    print(gammas[i], Wc)
