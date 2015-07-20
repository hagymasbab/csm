function [X,V,G,Z] = msmGenerate(N,randseed,A,B,sigma_x,sigma_v,g_shape,g_scale,z_shape,z_scale)
    % TODO test this
    setrandseed(randseed)
    Dv = size(A,2);
    Dx = size(A,1);
    k = size(B,2);
    
    G = gamrnd(g_shape,g_scale,[N k]);
    Z = gamrnd(z_shape,z_scale,[N 1]);
    
    V = mvnrnd(G * B',sigma_v*eye(Dv));
    
    Av = V  * A';
    X = mvnrnd(repmat(Z,1,Dx) .* Av,sigma_x*eye(Dx));
    
end