function [mu_post,C_post] = gsmPosteriorV(x,A,C,x_sigma,z_shape,z_scale,z_res)
    
    ATA = A' * A;
    sATA = ATA / x_sigma^2;    
    iC = stableInverse(C);
    sAx = A' * x / x_sigma^2;
    Dv = size(A,2);

    [z_post_dens,z_vals] = gsmPosteriorZ(x,A,C,x_sigma,z_shape,z_scale,z_res,true);
    
    component_Cs = zeros(z_res,Dv,Dv);
    component_mus = zeros(z_res,Dv);
    
    mu_post = zeros(Dv,1);
    C_post = zeros(Dv);
    
    for i=1:z_res
        Cc = stableInverse(iC + z_vals(i)^2 * sATA);
        component_Cs(i,:,:) = Cc;
        muc = z_vals(i) * Cc * sAx;
        component_mus(i,:) = muc;
        mu_post = mu_post + z_post_dens(i) * muc;
        C_post = C_post + z_post_dens(i) * Cc;
    end
    
    for i=1:z_res
        mu_diff = component_mus(i,:)' - mu_post;
        C_post = C_post + z_post_dens(i) * (mu_diff * mu_diff');
    end
end