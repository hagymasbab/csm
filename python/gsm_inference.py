import pickle
from numpy import zeros, identity, mean
from numpy.random import gamma, multivariate_normal
from pystan import StanModel

recompile = False
N = 2
d_x = 256
d_u = d_x
sigma_x = 0.7 * identity(d_x)
z_scale = 2
z_shape = 2
A = identity(d_x)
C = identity(d_u)

z_synth = gamma(z_shape, z_scale, N)
u_synth = multivariate_normal(zeros(d_u), C, N)
x_synth = zeros((N, d_x))
for i in range(N):
    act_u = u_synth[i, :].T
    act_mean = z_synth[i] * A.dot(act_u)
    x_synth[i, :] = multivariate_normal(act_mean, sigma_x)

gsm_dat = {
    'N': N,
    'd_x': d_x,
    'd_u': d_u,
    'sigma_x': sigma_x,
    'x': x_synth,
    'A': A,
    'C': C,
    'z_shape': z_shape,
    'z_scale': z_scale
}

if recompile:
    sm = StanModel(file='gsm_inference.stan')
    with open('gsm_inference.pkl', 'wb') as f:
        pickle.dump(sm, f)
else:
    sm = pickle.load(open('gsm_inference.pkl', 'rb'))

fit = sm.sampling(data=gsm_dat, iter=100, chains=8)
estimation = fit.extract(permuted=True)
z_est_mean = mean(estimation["z"], 0).T
print(z_est_mean)

print(z_synth)
