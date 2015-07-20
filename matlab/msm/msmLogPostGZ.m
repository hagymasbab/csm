function lp = msmLogPostGZ(g,z,x,A,B,sigma_x,sigma_v,g_shape,g_scale,z_shape,z_scale)
    Dv = size(A,2);
    AAT = A * A';
    Cxm = sigma_x * eye(Dv) + z^2 * sigma_v * AAT;
    
    lp = logGauss(x,z*A*B*g,Cxm) + logGamma(g,g_shape,g_scale) + logGamma(z,z_shape,z_scale);
    
end