import numpy as np
import helper as hp
import os
import multiprocessing as mp
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": True
    #"font.family": "Helvetica"
    })

sigma = np.array([[[1,0],[0,1]], [[0,1],[1,0]], [[0,-1j],[1j,0]], [[1,0],[0,-1]]])
#print(sigma)

def BBH3D(kz, ky, kx, gamma = 0.5, delta = 0, **kwargs):
    return np.kron(sigma[3], (np.kron(sigma[1], sigma[0])*(np.cos(kx) + gamma) - np.kron(sigma[2],sigma[3])*np.sin(kx))) - \
            np.kron(sigma[3], np.kron(sigma[2], sigma[2]*(np.cos(ky) + gamma) + sigma[1]*np.sin(ky))) + \
            np.kron((sigma[1]*(np.cos(kz) + gamma) - sigma[2]*np.sin(kz)), np.kron(sigma[0], sigma[0]))
            #np.kron(sigma[3], np.kron(sigma[3], sigma[0]))*delta
            #np.kron(sigma[3], np.kron(sigma[0] + sigma[3], sigma[0] + sigma[3]))*delta

def BBH3D_eff(kx, ky, kz, gamma = 0.5, delta = 0, Sigmax = 0, Sigmay = 0, Sigmaz = 0, **kwargs):
    return np.kron(sigma[3], (np.kron(sigma[1], sigma[0])*(np.cos(kx) + gamma + Sigmax) - np.kron(sigma[2],sigma[3])*np.sin(kx))) - \
            np.kron(sigma[3], np.kron(sigma[2], sigma[2]*(np.cos(ky) + gamma + Sigmay) + sigma[1]*np.sin(ky))) + \
            np.kron((sigma[1]*(np.cos(kz) + gamma + Sigmaz) - sigma[2]*np.sin(kz)), np.kron(sigma[0], sigma[0]))

#print( np.kron(sigma[0], np.kron(sigma[3] + sigma[0], sigma[0] + sigma[3])))
def test_ham():
    path = np.array([[0,0,0], [0,np.pi,0], [np.pi, np.pi, 0], [0,0,0], [np.pi, np.pi, np.pi], [0,np.pi,0]])
    n = 10
    bands = []
    args = {'gamma' : 2, 'delta' : 0.1}
    for i, point in enumerate(path[:-1]):
        for e in range(n):
            ks = point + e*(path[i+1]-point)/n
            bands.append(np.linalg.eigh(BBH3D(ks[0], ks[1], ks[2], **args))[0])
    bands = np.transpose(np.array(bands))

    fig, ax = plt.subplots()
    for band in bands:
        ax.plot(band)
    plt.show()


def wilson_loop(kx, ky, ham, args):
    vecs = []
    for kz in np.linspace(0, 2*np.pi, k_int)[:-1]:
        vecs.append(np.linalg.eigh(ham(kx, ky, kz, **args))[1][:,0:4])
    vecs = np.array(vecs)
    vecs2 = np.copy(vecs[1:])
    vecs2 = np.append(vecs2, [vecs[0]], axis = 0)
    loop = np.linalg.multi_dot(np.einsum('abc,abd->acd', np.conjugate(vecs), vecs2))
    return loop + args['delta']*np.kron(sigma[3] ,sigma[0])

def phases(vec):
    return np.real(-1j*np.log(vec))/(2*np.pi)

def test_wilson():
    args = {'gamma' : 2, 'delta' : 0.1}
    fig = plt.figure()
    ax = fig.add_subplot(projection = '3d')
    x,y,zs = [],[],[[],[],[],[]]
    for kx in np.linspace(0,2*np.pi,20):
        for ky in np.linspace(0,2*np.pi,20):
            x.append(kx)
            y.append(ky)
            eigs = np.linalg.eig(wilson_loop(kx, ky, BBH3D, args))[0]
            nus = np.sort(phases(eigs))
            for i,phase in enumerate(nus):
                zs[i].append(phase)

    for z in zs:
        ax.plot_trisurf(x,y,z)
    plt.show()

def wannier_bands(args):
    fig = plt.figure()
    ax = fig.add_subplot(projection = '3d')
    x,y,zs = [],[],[[],[],[],[]]
    for kx in np.linspace(-np.pi,np.pi,60):
        for ky in np.linspace(-np.pi,np.pi,60):
            x.append(kx)
            y.append(ky)
            eigs = np.linalg.eig(wilson_loop(kx, ky, BBH3D, args))[0]
            nus = np.sort(phases(eigs))
            for i,phase in enumerate(nus):
                zs[i].append(phase)

    for z in zs:
        ax.plot_trisurf(x,y,z)

    ax.set(xlabel = r'$k_y$', ylabel = r'$k_z$', zlabel = r'$\nu_x$')
    ax.set_zlabel(r'$\nu_x$', labelpad = 0)
    ax.tick_params(axis='x', which='major', pad=-3)
    ax.tick_params(axis='y', which='major', pad=-3)

    ax.view_init(10,-30)
    name = "WannierBandsBBH3D_gamma" + str(args['gamma'])
    #fig.set_size_inches(4, 3)
    fig.savefig(hp.plot_dir() + name + ".png", dpi = 300, bbox_inches='tight', pad_inches = 0.2)
    fig.savefig(hp.plot_dir() + name + ".pdf", bbox_inches = 'tight', pad_inches = 0.2)

def nested_wilson_loop(kx, ham, args):
    vecs = []
    for ky in np.linspace(0, 2*np.pi, k_int)[:-1]:
        eigvals, eigvecs = np.linalg.eig(wilson_loop(kx, ky, ham, args))
        sort = np.argsort(phases(eigvals))
        eigvecs = eigvecs[:,sort][:,:2]
        eigvecs2 = np.linalg.eigh(ham(kx, ky, 0, **args))[1][:,0:4]
        vecs.append(np.einsum('ac,ba', eigvecs, eigvecs2))
    vecs = np.array(vecs)
    vecs2 = np.copy(vecs[1:])
    vecs2 = np.append(vecs2, [vecs[0]], axis = 0)
    loop = np.linalg.multi_dot(np.einsum('abc,abd->acd', np.conjugate(vecs), vecs2))
    return loop, vecs[0]

def test_nested_wilson():
    args = {'gamma' : 0.5,'delta' : 0.01}
    nus = []
    x = []
    for kx in np.linspace(0,2*np.pi,50):
        x.append(kx)
        eigs = np.linalg.eig(nested_wilson_loop(kx, BBH3D, args)[0])[0]
        nus.append(np.sort(phases(eigs)))
    nus = np.transpose(np.array(nus))

    fig, ax = plt.subplots()
    for nu in nus:
        ax.plot(x, nu)
    plt.show()

def nested_nested_wilson_loop(ham, args):
    vecs = []
    for kx in np.linspace(0, 2*np.pi, k_int)[:-1]:
        loop, eigvecs2 = nested_wilson_loop(kx, ham, args)
        eigvals, eigvecs = np.linalg.eig(loop)
        sort = np.argsort(phases(eigvals))
        eigvecs = eigvecs[:,sort][:,:1]
        vecs.append(np.einsum('ac,ba', eigvecs, eigvecs2))
    vecs = np.array(vecs)
    vecs2 = np.copy(vecs[1:])
    vecs2 = np.append(vecs2, [vecs[0]], axis = 0)
    loop = np.linalg.multi_dot(np.einsum('abc,abd->acd', np.conjugate(vecs), vecs2))
    return phases(loop[0])[0]


def invariant():
    inv = []
    args = {'gamma' : 0, 'delta' : 0.00001}
    gammas = np.linspace(0, 2, 50)
    for gamma in gammas:
        args['gamma'] = gamma
        inv.append(np.abs(nested_nested_wilson_loop(BBH3D, args)))

    fig, ax = plt.subplots()
    #ax.set(xlabel = r'$\left|\frac{\gamma}{\lambda}\right|$', ylabel = r'$q_{xy}$', fontsize = 12)
    ax.set_xlabel(r'$\left|\frac{\gamma}{\lambda}\right|$', fontsize = 14)
    ax.set_ylabel(r'$p_z$', fontsize = 14)
    ax.plot(gammas, inv, ls = '', marker = '.', markersize = 4)

    name = "InvBBH3D"
    fig.set_size_inches(4, 3)
    fig.savefig(hp.plot_dir() + name + ".png", dpi = 300, bbox_inches = 'tight', pad_inches = 0)
    fig.savefig(hp.plot_dir() + name + ".pdf", bbox_inches = 'tight', pad_inches = 0)


k_int = 20

#test_ham()
#test_wilson()
#test_nested_wilson()

#args = {'gamma' : 1.1, 'delta' : 0.00001, 'Sigmax' : 0.2, 'Sigmay' : -2, 'Sigmaz' : -2}
#print(np.abs(nested_nested_wilson_loop(BBH3D_eff, args)))

args = {'gamma' : 0.5, 'delta' : 0.0}
wannier_bands(args)
args['gamma'] = 1
wannier_bands(args)
args['gamma'] = 2
wannier_bands(args)

#invariant()
