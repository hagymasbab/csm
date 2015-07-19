function grad = msmLogPostGZGrad(x,z,g,A,B,sigma_x,sigma_v,g_shape,g_scale,z_shape,z_scale)
    Dv = size(A,2);
    AAT = A * A';
    Cxm = sigma_x * eye(Dv) + z^2 * sigma_v * AAT;
    iCxm = stableInverse(Cxm);
    zAB = z * A * B;
    y = zAB * g;
    
    g_grad = zAB' * iCxm * (x - y) + logGammaDerivative(g,g_shape,g_scale);
    
    xx = x * x';
    xy = x * y';
    yy = y * y';
    z_mat = iCxm * ( z * sigma_v * AAT * (xx + z * xy) - sigma_x * (xy + z * yy) ) * iCxm;
    z_grad = -2 * trace(z_mat) + logGammaDerivative(z,z_shape,z_scale);
    
    grad = [g_grad;z_grad];
    
end