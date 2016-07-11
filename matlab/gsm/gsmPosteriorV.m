function [mu_post,C_post,z_mean,z_var,z_post_dens,component_mus,component_Cs] = gsmPosteriorV(x,A,C,x_sigma,z_shape,z_scale,z_res)
    
    ATA = A' * A;
    sATA = ATA / x_sigma^2;    
    iC = stableInverse(C);
    sAx = A' * x / x_sigma^2;
    Dv = size(A,2);
    
    fprintf('z')
    [z_post_dens,z_vals] = gsmPosteriorZ(x,A,C,x_sigma,z_shape,z_scale,z_res,true);
    
%     size(z_post_dens)
%     size(z_vals)
    
    z_mean = z_post_dens' * z_vals;
    z_var = z_post_dens' * (z_vals - z_mean).^2;
    
    %component_Cs = zeros(z_res,Dv,Dv);
    component_Cs = cell(1,z_res);
    %component_mus = zeros(z_res,Dv);
    component_mus = cell(1,z_res);
    
    mu_post = zeros(Dv,1);
    C_post = zeros(Dv);
    
    fprintf('int')
    
    for i=1:z_res
        Cc = stableInverse(iC + z_vals(i)^2 * sATA);
        %component_Cs(i,:,:) = Cc;
        component_Cs{i} = Cc;
        muc = z_vals(i) * Cc * sAx;
        %component_mus(i,:) = muc;
        component_mus{i} = muc;
        mu_post = mu_post + z_post_dens(i) * muc;
        C_post = C_post + z_post_dens(i) * Cc;
    end
    
    for i=1:z_res
        %mu_diff = component_mus(i,:)' - mu_post;
        mu_diff = component_mus{i} - mu_post;
        C_post = C_post + z_post_dens(i) * (mu_diff * mu_diff');
    end
end