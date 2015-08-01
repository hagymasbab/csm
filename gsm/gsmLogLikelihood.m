function loglike = gsmLogLikelihood(X,C,A,x_sigma,z_shape,z_scale,z_res)
    
    Du = size(C,1);
    Dx = size(X,2);
    if size(A,1) ~= Dx
        error('Data is of dimension %d and filterbank output dimension is %d',Dx,size(A,1));
    end
    if size(A,2) ~= Du
        error('Latent covariance dimension %d and filterbank intput dimension is %d',Du,size(A,2));
    end
    
    iC = inv(C);        
    noiseCov = eye(Dx) * x_sigma^2;
    x0 = zeros(Dx,1);
    ACAT = A * C * A';        
    
    loglike = 0;    
    for t=1:size(X,1)
        printCounter(t,'maxVal',size(X,1),'stringVal','Observation');
        xt = X(t,:)';
        
        % get a MAP estimate of z
        [z_map,z_min,z_max] = get_pz_x_max(xt, noiseCov, ACAT, x0, z_shape, z_scale);
        z_vals = linspace(z_min,z_max,z_res);
        z_square = z_vals.^2;
        
        act_like = 0;
        for i = 1:z_res
            like_cov = noiseCov + z_square(i) * ACAT;
            likelihood = mvnpdf(xt,zeros(size(xt)),like_cov);
            prior = gampdf(z_vals(i),z_shape,z_scale);
            act_like = act_like + likelihood * prior;
        end
        
        loglike = loglike + log(act_like);
        
    end
end