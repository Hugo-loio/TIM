import numpy as np
import helper as hp
import os
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": True
    #"font.family": "Helvetica"
    })

sigma = np.array([[[1,0],[0,1]], [[0,1],[1,0]], [[0,-1j],[1j,0]], [[1,0],[0,-1]]])
#print(sigma)

def BBH2D(ky, kx, gamma = 0.5, delta = 0, **kwargs):
    return np.kron(sigma[1], sigma[0])*(np.cos(kx) + gamma) \
        - np.kron(sigma[2], sigma[2])*(np.cos(ky) + gamma) \
        - np.kron(sigma[2], (sigma[3]*np.sin(kx) + sigma[1]*np.sin(ky))) \
        + delta*np.kron(sigma[3] ,sigma[0])

def wilson_loop(kx, ham, args):
    vecs = []
    for ky in np.linspace(0, 2*np.pi, k_int)[:-1]:
        vecs.append(np.linalg.eigh(ham(kx, ky, **args))[1][:,0:2])
    vecs = np.array(vecs)
    vecs2 = np.copy(vecs[1:])
    vecs2 = np.append(vecs2, [vecs[0]], axis = 0)
    loop = np.linalg.multi_dot(np.einsum('abc,abd->acd', np.conjugate(vecs), vecs2))
    return loop 

def phases(vec):
    return np.real(-1j*np.log(vec))/(2*np.pi)

def nested_wilson_loop(ham, args):
    vecs = []
    for kx in np.linspace(0, 2*np.pi, k_int)[:-1]:
        eigvals, eigvecs = np.linalg.eig(wilson_loop(kx, ham, args))
        sort = np.argsort(phases(eigvals))
        eigvecs = eigvecs[:,sort][:,:1]
        eigvecs2 = np.linalg.eigh(ham(kx, 0, **args))[1][:,0:2]
        vecs.append(np.einsum('ac,ba', eigvecs, eigvecs2))
    vecs = np.array(vecs)
    vecs2 = np.copy(vecs[1:])
    vecs2 = np.append(vecs2, [vecs[0]], axis = 0)
    loop = np.linalg.multi_dot(np.einsum('abc,abd->acd', np.conjugate(vecs), vecs2))
    return phases(loop[0])[0]

def wannier_bands(args):
    nus = []
    x = []
    for kx in np.linspace(-np.pi,np.pi,401):
        x.append(kx)
        eigs = np.linalg.eig(wilson_loop(kx, BBH2D, args))[0]
        nus.append(np.sort(phases(eigs)))
    nus = np.transpose(np.array(nus))

    fig, ax = plt.subplots()
    ax.set(xlabel = r'$k_y$', ylabel = r'$\nu$')
    for nu in nus:
        ax.plot(x, nu)

    name = "WannierBandsBBH2D_gamma" + str(args['gamma'])
    fig.set_size_inches(3.4, 2.5)
    fig.savefig(hp.plot_dir() + name + ".png", dpi = 300, bbox_inches = 'tight', pad_inches = 0)
    fig.savefig(hp.plot_dir() + name + ".pdf", bbox_inches = 'tight', pad_inches = 0)
    #plt.show()

def invariant():
    inv = []
    args = {'gamma' : 0}
    gammas = np.linspace(0, 2, 400)
    for gamma in gammas:
        args['gamma'] = gamma
        inv.append(np.abs(nested_wilson_loop(BBH2D, args)))

    fig, ax = plt.subplots()
    #ax.set(xlabel = r'$\left|\frac{\gamma}{\lambda}\right|$', ylabel = r'$q_{xy}$', fontsize = 12)
    ax.set_xlabel(r'$\left|\frac{\gamma}{\lambda}\right|$', fontsize = 14)
    ax.set_ylabel(r'$p_y^{\nu_x^-}$', fontsize = 14)
    ax.plot(gammas, inv, ls = '--', marker = '.', markersize = 3)

    name = "InvBBH2D"
    fig.set_size_inches(4, 3)
    fig.savefig(hp.plot_dir() + name + ".png", dpi = 300, bbox_inches = 'tight', pad_inches = 0)
    fig.savefig(hp.plot_dir() + name + ".pdf", bbox_inches = 'tight', pad_inches = 0)
    #plt.show()

k_int = 30
                  
args = {'gamma' : 0.5, 'delta' : 0.0}
print(np.abs(nested_wilson_loop(BBH2D, args)))

#wannier_bands(args)
args['gamma'] = 1
#wannier_bands(args)
args['gamma'] = 2
#wannier_bands(args)

invariant()
