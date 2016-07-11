function [z_post_dens,z_vals] = gsmPosteriorZ(x,A,C,x_sigma,z_shape,z_scale,z_res,scinot)
    ACAT = A * C * A';
    noiseCov = eye(size(x,1)) * x_sigma^2;
    
    % get a MAP estimate of z
    %[~,z_min,z_max] = get_pz_x_max(x, noiseCov, ACAT, zeros(size(x)), z_shape, z_scale);
    z_min = 0.1;z_max = 5;
    z_vals = linspace(z_min,z_max,z_res)';    
    
    fprintf('zmapd')
    
    if ~scinot
        z_post_numerators = zeros(z_res,1);    
        for i = 1:z_res
            like_cov = noiseCov + z_vals(i)^2 * ACAT;
            likelihood = stableMvnpdf(x,zeros(size(x)),like_cov);
            prior = gampdf(z_vals(i),z_shape,z_scale);
            z_post_numerators(i,1) = likelihood * prior;
        end
        z_post_denominator = sum(z_post_numerators);
        z_post_dens = z_post_numerators / z_post_denominator;    
    else
        z_post_numerators_expo = zeros(z_res,1);
        z_post_numerators_coeff = zeros(z_res,1);
        z_post_denominator_coeff = 0;
        z_post_denominator_expo = 0;
        for i = 1:z_res
            like_cov = noiseCov + z_vals(i)^2 * ACAT;
            [like_coeff,like_expo] = stableMvnpdf(x,zeros(size(x)),like_cov,'scientific',true);
            [prior_coeff,prior_expo] = sciNot(gampdf(z_vals(i),z_shape,z_scale),false);
            [z_post_numerators_coeff(i),z_post_numerators_expo(i)] = prodSciNot([like_coeff prior_coeff],[like_expo prior_expo]);
            [z_post_denominator_coeff,z_post_denominator_expo] = sumSciNot(z_post_denominator_coeff,z_post_denominator_expo,z_post_numerators_coeff(i),z_post_numerators_expo(i));        
        end
        z_post_dens = zeros(z_res,1);
        for i = 1:z_res
            [z_post_coeff,z_post_expo] = prodSciNot([z_post_numerators_coeff(i) 1/z_post_denominator_coeff],[z_post_numerators_expo(i) -z_post_denominator_expo]);
            z_post_dens(i) = z_post_coeff * 10^z_post_expo;
        end
    end
end