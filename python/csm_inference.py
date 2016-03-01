import pickle
from numpy import zeros, ones, identity, mean, minimum
from numpy.random import gamma, dirichlet, multivariate_normal
from pystan import StanModel

recompile = False
model_type = 2
N = 1
k = 2
d_x = 16
d_u = d_x
sigma_x = 0.7 * identity(d_x)
z_scale = 2
z_shape = 2
g_scale = 0.1
g_shape = 1
A = identity(d_x)
Var_u = identity(d_x)
templates = [[0, 1, 2, 3, 4, 5, 6, 7], [8, 9, 10, 11, 12, 13, 14, 15]]
C = []
for i in range(k):
    if model_type == 1:
        act_C = identity(d_u)
    elif model_type == 2:
        act_C = zeros((d_u, d_u))
    for j in range(len(templates[i])):
        idx1 = templates[i][j]
        for h in range(j+1, len(templates[i])):
            idx2 = templates[i][h]
            act_C[idx1, idx2] = 0.9
            act_C[idx2, idx1] = act_C[idx1, idx2]
    C.append(act_C)

z_synth = gamma(z_shape, z_scale, N)
g_synth = gamma(g_shape, g_scale, (N, k))
# g_synth = dirichlet((1, 1), N)
# u_synth = multivariate_normal(zeros(d_u), C, N)
x_synth = zeros((N, d_x))
for i in range(N):
    act_C = zeros((d_u, d_u))
    for j in range(k):
        act_C += g_synth[i, j] * C[j]
    if model_type == 2:
        act_C = minimum(act_C, ones((d_u, d_u))) + identity(d_u)
        act_C = Var_u * act_C * Var_u
    act_u = multivariate_normal(zeros(d_u), act_C)
    act_mean = z_synth[i] * A.dot(act_u)
    x_synth[i, :] = multivariate_normal(act_mean, sigma_x)

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

if recompile:
    sm = StanModel(file=fname + '_inference.stan')
    with open(fname + 'inference.pkl', 'wb') as f:
        pickle.dump(sm, f)
else:
    sm = pickle.load(open(fname + '_inference.pkl', 'rb'))

fit = sm.sampling(data=gsm_dat, iter=400, chains=2)
estimation = fit.extract(permuted=True)

g_est_mean = mean(estimation["g"], 0)
print('g est', g_est_mean)
print('g true', g_synth)

z_est_mean = mean(estimation["z"], 0).T
print('z est', z_est_mean)
print('z true', z_synth)
