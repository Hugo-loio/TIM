import numpy as np
import os
import multiprocessing as mp

sigma = np.array([[[1,0],[0,1]], [[0,1],[1,0]], [[0,-1j],[1j,0]], [[1,0],[0,-1]]])
#print(sigma)

Gamma0 = np.kron(sigma[3],(np.kron(sigma[1], sigma[0]) - np.kron(sigma[2], sigma[2]))) + np.kron(sigma[1], np.kron(sigma[0], sigma[0]))
#print(np.real(Gamma0))

alpha = [[1,3],[1,4],[1,5],[2,3],[2,4],[2,6],[3,7],[4,8],[5,7],[5,8],[6,7],[6,8]]

Us = []
for i in range(12):
    Us.append(np.zeros(Gamma0.shape))
    ind = np.array(alpha[i])-1
    Us[-1][ind[0],ind[1]] = 1
    Us[-1][ind[1],ind[0]] = 1

#print(U[2])
def H0(kx, ky, kz, gamma):
    return np.kron(sigma[3], (np.kron(sigma[1], sigma[0])*(np.cos(kx) + gamma) - np.kron(sigma[2],sigma[3])*np.sin(kx))) - \
            np.kron(sigma[3], np.kron(sigma[2], sigma[2]*(np.cos(ky) + gamma) + sigma[1]*np.sin(kx))) + \
            np.kron((sigma[1]*(np.cos(kz) + gamma) - sigma[2]*np.sin(kz)), np.kron(sigma[0], sigma[0]))

def aux_disc(gamma, ndisc): 
    ks = np.linspace(0,2*np.pi,ndisc)
    return np.array([(E+1j*1E-9)*np.identity(8) - H0(kx,ky,kz,gamma) for kx in ks for ky in ks for kz in ks])

def G(kx, ky, kz, gamma, E, Sigma):
    return np.linalg.inv((E+1j*1E-9)*np.identity(8)-H0(kx,ky,kz,gamma)-Sigma)

def update_Sigma(Sigma, gamma, E, W, aux):
    #int_res = [np.linalg.inv(M-Sigma) for M in aux]
    int_res = np.average(np.linalg.inv(aux-Sigma), axis = 0)
    res = np.zeros((8,8), dtype='complex128')
    for U in Us:
        res += np.linalg.multi_dot([U, int_res, U])
    res *= W**2/12
    return res

def renorm_sigmas(Sigma0, E, gamma, W, aux):
    Sigma1 = Sigma0.copy()
    max_order = 100
    for i in range(max_order):
        Sigma2 = update_Sigma(Sigma1, gamma, E, W, aux)
        norm = np.linalg.norm(Sigma2-Sigma1)
        Sigma1 = Sigma2
        if(norm < 1E-5):
            break
    if(i == max_order-1):
        print("ERROR: renormalization did not converge",gamma, W)
    if(W == W_max):
        print("Did gamma", gamma)
    return np.real(np.array([Sigma1[0,2], Sigma1[0,3], Sigma1[0,4]])) + gamma

E = 0
ndisc = 10
Sigma0 = np.zeros((8,8))
gamma_min = 1
gamma_max = 1.3
W_min = 0
W_max = 5
ngrid = 200
version = 1

gammas = np.linspace(gamma_min, gamma_max, ngrid)
Ws = np.linspace(W_min, W_max, ngrid)

print("Generating aux matrices...")
pool = mp.Pool(mp.cpu_count(), initargs= (ndisc))
auxs = [pool.apply_async(aux_disc, args= (gamma, ndisc)) for gamma in gammas]
pool.close()
pool.join()
auxs = [aux.get() for aux in auxs]

print("Renormalizing...")

pool = mp.Pool(mp.cpu_count())
res = [[pool.apply_async(renorm_sigmas, args=(Sigma0, E, gamma, W, auxs[i])) for W in Ws] for i,gamma in enumerate(gammas)]
#res = [pool.apply(renorm_sigmas, args=(Sigma0, E, gammas[0], W, auxs[0])) for W in Ws]
pool.close()
pool.join()
res = np.array([[r.get() for r in re] for re in res])

file_name = "SCBA_ndisc" + str(ndisc) + "_v" + str(version) 
path_name = os.path.dirname(os.path.realpath(__file__)) + "/data/" + file_name
np.save(path_name, res)

file_name = "SCBA_ndisc" + str(ndisc) + "_v" + str(version) + "_params"
path_name = os.path.dirname(os.path.realpath(__file__)) + "/data/" + file_name
np.save(path_name, np.array([Ws,gammas]))

#res = np.amax(res, axis = 2)
