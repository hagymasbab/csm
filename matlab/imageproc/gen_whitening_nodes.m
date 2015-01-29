function wnode = gen_whitening_nodes(x)

x = reshape(x, size(x,1), size(x,2)*size(x,3));

% calculate the sample mean and covariance matrix for the data
% x_bar = mean(x);
sigma = cov(x);
% Get the eigenvectors and corresponding eigenvalues for this data.
[U, lambda] = eig(sigma);
% Whiten the data by applying the following linear transform. and plot the
% result.

wnode = U*lambda^(-.5)*U';
% wnode = lambda^(-.5)*U;

