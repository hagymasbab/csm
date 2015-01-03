function [samples, epsilon] = nuts_da(f, M, Madapt, theta0, delta)
% function [samples, epsilon] = nuts_da(f, M, Madapt, theta0, delta)
% 
% Implements the No-U-Turn Sampler (NUTS), specifically, algorithm 6
% from the NUTS paper (Hoffman & Gelman, 2011). Runs Madapt steps of
% burn-in, during which it adapts the step size parameter epsilon, then
% starts generating samples to return.
% 
% epsilon is a step size parameter.
% f(theta) should be a function that returns the log probability its
% gradient evaluated at theta. I.e., you should be able to call 
% [logp grad] = f(theta).
% M is the number of samples to generate.
% Madapt is the number of steps of burn-in/how long to run the dual
% averaging algorithm to fit the step size epsilon. 
% theta0 is a 1-by-D vector with the desired initial setting of the parameters.
% delta should be between 0 and 1, and is a target HMC acceptance
% probability. Defaults to 0.6 if unspecified.
% 
% The returned variable "samples" is an M-by-D matrix of post-burn-in samples 
% generated by NUTS.
% The returned variable "epsilon" is the step size that was fit using dual
% averaging.

% Copyright (c) 2011, Matthew D. Hoffman
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 
% Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

global nfevals;
nfevals = 0;

if (nargin < 5)
    delta = 0.6;
end

assert(size(theta0, 1) == 1);

D = length(theta0);
samples = zeros(M+Madapt, D);

[logp grad] = f(theta0);
samples(1, :) = theta0;

% Choose a reasonable first epsilon by a simple heuristic.
epsilon = find_reasonable_epsilon(theta0, grad, logp, f);

% Parameters to the dual averaging algorithm.
gamma = 0.05;
t0 = 10;
kappa = 0.75;
mu = log(10*epsilon);
% Initialize dual averaging algorithm.
epsilonbar = 1;
Hbar = 0;

for m = 2:M+Madapt,
    % Resample momenta.
    r0 = randn(1, D);
    % Joint log-probability of theta and momenta r.
    joint = logp - 0.5 * (r0 * r0');
    % Resample u ~ uniform([0, exp(joint)]).
    % Equivalent to (log(u) - joint) ~ exponential(1).
    logu = joint - exprnd(1);
    % Initialize tree.
    thetaminus = samples(m-1, :);
    thetaplus = samples(m-1, :);
    rminus = r0;
    rplus = r0;
    gradminus = grad;
    gradplus = grad;
    % Initial height j = 0.
    j = 0;
    % If all else fails, the next sample is the previous sample.
    samples(m, :) = samples(m-1, :);
    % Initially the only valid point is the initial point.
    n = 1;

    % Main loop---keep going until the criterion s == 0.
    s = 1;
    while (s == 1)
        % Choose a direction. -1=backwards, 1=forwards.
        v = 2*(rand() < 0.5)-1;
        % Double the size of the tree.
        if (v == -1)
            [thetaminus, rminus, gradminus, ~, ~, ~, thetaprime, gradprime, logpprime, nprime, sprime, alpha, nalpha] = ...
                build_tree(thetaminus, rminus, gradminus, logu, v, j, epsilon, f, joint);
        else
            [~, ~, ~, thetaplus, rplus, gradplus, thetaprime, gradprime, logpprime, nprime, sprime, alpha, nalpha] = ...
                build_tree(thetaplus, rplus, gradplus, logu, v, j, epsilon, f, joint);
        end
        % Use Metropolis-Hastings to decide whether or not to move to a
        % point from the half-tree we just generated.
        if ((sprime == 1) && (rand() < nprime/n))
            samples(m, :) = thetaprime;
            logp = logpprime;
            grad = gradprime;
        end
        % Update number of valid points we've seen.
        n = n + nprime;
        % Decide if it's time to stop.
        s = sprime && stop_criterion(thetaminus, thetaplus, rminus, rplus);
        % Increment depth.
        j = j + 1;
    end
    
    % Do adaptation of epsilon if we're still doing burn-in.
    eta = 1 / (m - 1 + t0);
    Hbar = (1 - eta) * Hbar + eta * (delta - alpha / nalpha);
    if (m <= Madapt)
        epsilon = exp(mu - sqrt(m-1)/gamma * Hbar);
        eta = (m-1)^-kappa;
        epsilonbar = exp((1 - eta) * log(epsilonbar) + eta * log(epsilon));
    else
        epsilon = epsilonbar;
    end
end
samples = samples(Madapt+1:end, :);
fprintf('Took %d gradient evaluations.\n', nfevals);
end

function [thetaprime, rprime, gradprime, logpprime] = leapfrog(theta, r, grad, epsilon, f)
rprime = r + 0.5 * epsilon * grad;
thetaprime = theta + epsilon * rprime;
[logpprime, gradprime] = f(thetaprime);
rprime = rprime + 0.5 * epsilon * gradprime;
global nfevals;
nfevals = nfevals + 1;
end

function criterion = stop_criterion(thetaminus, thetaplus, rminus, rplus)
thetavec = thetaplus - thetaminus;
criterion = (thetavec * rminus' >= 0) && (thetavec * rplus' >= 0);
end

% The main recursion.
function [thetaminus, rminus, gradminus, thetaplus, rplus, gradplus, thetaprime, gradprime, logpprime, nprime, sprime, alphaprime, nalphaprime] = ...
                build_tree(theta, r, grad, logu, v, j, epsilon, f, joint0)
if (j == 0)
    % Base case: Take a single leapfrog step in the direction v.
    [thetaprime, rprime, gradprime, logpprime] = leapfrog(theta, r, grad, v*epsilon, f);
    joint = logpprime - 0.5 * (rprime * rprime');
    % Is the new point in the slice?
    nprime = logu < joint;
    % Is the simulation wildly inaccurate?
    sprime = logu - 1000 < joint;
    % Set the return values---minus=plus for all things here, since the
    % "tree" is of depth 0.
    thetaminus = thetaprime;
    thetaplus = thetaprime;
    rminus = rprime;
    rplus = rprime;
    gradminus = gradprime;
    gradplus = gradprime;
    % Compute the acceptance probability.
    alphaprime = min(1, exp(logpprime - 0.5 * (rprime * rprime') - joint0));
    nalphaprime = 1;
else
    % Recursion: Implicitly build the height j-1 left and right subtrees.
    [thetaminus, rminus, gradminus, thetaplus, rplus, gradplus, thetaprime, gradprime, logpprime, nprime, sprime, alphaprime, nalphaprime] = ...
                build_tree(theta, r, grad, logu, v, j-1, epsilon, f, joint0);
    % No need to keep going if the stopping criteria were met in the first
    % subtree.
    if (sprime == 1)
        if (v == -1)
            [thetaminus, rminus, gradminus, ~, ~, ~, thetaprime2, gradprime2, logpprime2, nprime2, sprime2, alphaprime2, nalphaprime2] = ...
                build_tree(thetaminus, rminus, gradminus, logu, v, j-1, epsilon, f, joint0);
        else
            [~, ~, ~, thetaplus, rplus, gradplus, thetaprime2, gradprime2, logpprime2, nprime2, sprime2, alphaprime2, nalphaprime2] = ...
                build_tree(thetaplus, rplus, gradplus, logu, v, j-1, epsilon, f, joint0);
        end
        % Choose which subtree to propagate a sample up from.
        if (rand() < nprime2 / (nprime + nprime2))
            thetaprime = thetaprime2;
            gradprime = gradprime2;
            logpprime = logpprime2;
        end
        % Update the number of valid points.
        nprime = nprime + nprime2;
        % Update the stopping criterion.
        sprime = sprime && sprime2 && stop_criterion(thetaminus, thetaplus, rminus, rplus);
        % Update the acceptance probability statistics.
        alphaprime = alphaprime + alphaprime2;
        nalphaprime = nalphaprime + nalphaprime2;
    end
end
end

function epsilon = find_reasonable_epsilon(theta0, grad0, logp0, f)
epsilon = 1;
r0 = randn(1, length(theta0));
% Figure out what direction we should be moving epsilon.
[~, rprime, ~, logpprime] = leapfrog(theta0, r0, grad0, epsilon, f);
acceptprob = exp(logpprime - logp0 - 0.5 * (rprime * rprime' - r0 * r0'));
a = 2 * (acceptprob > 0.5) - 1;
% Keep moving epsilon in that direction until acceptprob crosses 0.5.
while (acceptprob^a > 2^(-a))
    epsilon = epsilon * 2^a;
    [~, rprime, ~, logpprime] = leapfrog(theta0, r0, grad0, epsilon, f);
    acceptprob = exp(logpprime - logp0 - 0.5 * (rprime * rprime' - r0 * r0'));
end
end









