function [z_post_dens,z_vals] = gsmPosteriorZ(x,A,C,x_sigma,z_shape,z_scale,z_res)
    ACAT = A * C * A';
    noiseCov = eye(size(x,1)) * x_sigma^2;
    
    % get a MAP estimate of z
    [~,z_min,z_max] = get_pz_x_max(x, noiseCov, ACAT, zeros(size(x)), z_shape, z_scale);
    z_vals = linspace(z_min,z_max,z_res);
    
    z_post_numerators = zeros(z_res,1);
    for i = 1:z_res
        like_cov = noiseCov + z_vals(i)^2 * ACAT;
        likelihood = stableMvnpdf(x,zeros(size(x)),like_cov);
        prior = gampdf(z_vals(i),z_shape,z_scale);
        z_post_numerators(i,1) = likelihood * prior;
    end
    z_post_denominator = sum(z_post_numerators);
    z_post_dens = z_post_numerators / z_post_denominator;
end