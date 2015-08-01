function [mu_post,C_post] = gsmPosteriorV(x,A,C,x_sigma,z_shape,z_scale,z_res)
    z_post_dens = gsmPosteriorZ(x,A,C,x_sigma,z_shape,z_scale,z_res);
end