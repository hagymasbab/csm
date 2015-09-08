function samples = sampleFromGaussianMixture(mixing,means,covs,N)
    Dx = length(means{1});
    samples = zeros(N,Dx);
    for i=1:N
        act_z = find(mnrnd(1,mixing(:)'));
        samples(i,:) = mvnrnd(means{act_z}',covs{act_z});
    end
end