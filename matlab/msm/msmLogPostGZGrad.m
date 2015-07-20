function grad = msmLogPostGZGrad(g,z,x,A,B,sigma_x,sigma_v,g_shape,g_scale,z_shape,z_scale)
    Dv = size(A,2);
    AAT = A * A';
    
    g = real(g);
    z = real(z);
    
     if ~isreal(g)
         g
        error('complex g!!');
    end
    
    if ~isreal(z)
        z
        error('complex z!!');
    end
    
    Cxm = sigma_x * eye(Dv) + z^2 * sigma_v * AAT;       
    iCxm = stableInverse(Cxm);
    
    zAB = z * A * B;
    y = zAB * g;
    
    g_grad = zAB' * iCxm * (x - y) + logGammaDerivative(g,g_shape,g_scale);
    
    if ~isreal(g_grad)
        error('complex g_grad!!');
    end
    
    xx = x * x';
    xy = x * y';
    yy = y * y';
    z_mat = iCxm * ( z * sigma_v * (xx - z * xy) * AAT + sigma_x * (xy - z * yy) ) * iCxm;
    z_grad = trace(z_mat) - z * sigma_v * trace(iCxm * AAT) + logGammaDerivative(z,z_shape,z_scale);
    
    if ~isreal(z_grad)
        error('complex z_grad!!');
    end
    
    grad = [g_grad;z_grad];
end