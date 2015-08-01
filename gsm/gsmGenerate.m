function X = gsmGenerate(N,C,A,x_sigma,z_shape,z_scale)
    Du = size(C,1);
    Dx = size(A,1);

    Z = gamrnd(z_shape,z_scale,N,1);
    U = mvnrnd(zeros(N,Du),C);
    
    x_mean = repmat(Z,1,Dx) .* (A * U')';
    
    X = mvnrnd(x_mean,x_sigma^2 * eye(Dx));
end