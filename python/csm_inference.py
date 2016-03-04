import pickle
from scipy.io import loadmat
import numpy as np
import numpy.random as rnd
from pystan import StanModel
import matplotlib.pyplot as pl

recompile = False
model_type = 2

N = 1
k = 2
d_x = 16
d_u = d_x
sigma_x = 0.1 * np.identity(d_x)
z_scale = 2
z_shape = 2
g_scale = 1
g_shape = 1
Var_u = np.identity(d_x)

A = np.identity(d_x)
i, j = np.indices(A.shape)
A[i == j+1] = 0.5
A[i == j-1] = 0.5
# matDict = loadmat('../matlab/bin/filters_OF_64.mat')
# A = matDict['A']
# A = A + 1 * np.identity(d_x)
# pl.imshow(A, interpolation='nearest')
# pl.show()

templates = [[0, 1, 2, 3, 4, 5, 6, 7], [8, 9, 10, 11, 12, 13, 14, 15]]
if model_type == 3:
    C = np.zeros((d_u, k))
else:
    C = []
for i in range(k):
    if model_type == 3:
        for j in range(len(templates[i])):
            C[templates[i][j], i] = 0.9
    else:
        if model_type == 1:
            act_C = np.identity(d_u)
        elif model_type == 2:
            act_C = np.zeros((d_u, d_u))
        for j in range(len(templates[i])):
            idx1 = templates[i][j]
            for h in range(j+1, len(templates[i])):
                idx2 = templates[i][h]
                act_C[idx1, idx2] = 0.9
                act_C[idx2, idx1] = act_C[idx1, idx2]
        C.append(act_C)

z_synth = rnd.gamma(z_shape, z_scale, N)

if model_type == 2:
    g_synth = rnd.beta(g_shape, g_scale, (N, k))
else:
    g_synth = rnd.gamma(g_shape, g_scale, (N, k))
# g_synth = dirichlet((1, 1), N)
u_synth = np.zeros((N, d_u))
x_synth = np.zeros((N, d_x))
for i in range(N):
    if model_type == 3:
        act_g = g_synth[i, :].T
        mu_u = np.dot(C, act_g)
        # mu_u = C.T.dot(act_g)
        act_u = rnd.multivariate_normal(mu_u, Var_u)
    else:
        act_C = np.zeros((d_u, d_u))
        for j in range(k):
            act_C += g_synth[i, j] * C[j]
        if model_type == 2:
            act_C = np.minimum(act_C, np.ones((d_u, d_u))) + np.identity(d_u)
            act_C = np.dot(np.dot(Var_u, act_C), Var_u)
            # pl.imshow(act_C, interpolation='nearest')
            # pl.show()
        act_u = rnd.multivariate_normal(np.zeros(d_u), act_C)
    u_synth[i, :] = act_u.T
    act_mean = z_synth[i] * A.dot(act_u)
    x_synth[i, :] = rnd.multivariate_normal(act_mean, sigma_x)

gsm_dat = {
    'N': N,
    'k': k,
    'd_x': d_x,
    'd_u': d_u,
    'sigma_x': sigma_x,
    'x': x_synth,
    'A': A,
    'C': C,
    'z_shape': z_shape,
    'z_scale': z_scale,
    'g_shape': g_shape,
    'g_scale': g_scale,
    'Var_u': Var_u
}

if model_type == 1:
    fname = "csm"
elif model_type == 2:
    fname = "csm2"
elif model_type == 3:
    fname = "msm"
if recompile:
    sm = StanModel(file=fname + '_inference.stan')
    with open(fname + '_inference.pkl', 'wb') as f:
        pickle.dump(sm, f)
else:
    sm = pickle.load(open(fname + '_inference.pkl', 'rb'))

fit = sm.sampling(data=gsm_dat, iter=2000, chains=2)
estimation = fit.extract(permuted=True)

g_est_mean = np.mean(estimation["g"], 0)
print('g est', g_est_mean)
print('g true', g_synth)

z_est_mean = np.mean(estimation["z"], 0).T
print('z est', z_est_mean)
print('z true', z_synth)

pl.subplot(221)
pl.hist(estimation['z'], bins=40)
pl.plot([z_synth[0], z_synth[0]], [0, pl.gca().get_ylim()[1]], color='r', linestyle='-', linewidth=2)

pl.subplot(222)
pl.hist2d(estimation['g'][:, 0, 0], estimation['g'][:, 0, 1], bins=40)
maxlim = max(pl.gca().get_xlim()[1], pl.gca().get_ylim()[1])
pl.gca().set_xlim([0, maxlim])
pl.gca().set_ylim([0, maxlim])
pl.scatter(g_synth[0, 0], g_synth[0, 1], marker='x', c='w', s=100, linewidth=2)
pl.plot([0, maxlim], [0, maxlim], color='w', linestyle='-', linewidth=2)

pl.subplot(223)
imdim = np.sqrt(d_u)
pl.imshow(u_synth[0, :].reshape((imdim, imdim)), interpolation="nearest")

pl.subplot(224)
u_est_mean = np.mean(estimation["u"], 0)
pl.imshow(u_est_mean.reshape((imdim, imdim)), interpolation="nearest")

pl.show()
