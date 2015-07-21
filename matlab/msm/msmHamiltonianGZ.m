function [gsamp,zsamp] = msmHamiltonianGZ(x,L,burnin,randseed,A,B,sigma_x,sigma_v,g_shape,g_scale,z_shape,z_scale,target_acceptance)
    
    setrandseed(randseed);
    Dv = size(A,2);
    k = size(B,2);
    x = reshape(x,Dv,1);    
    
    % set initial values for g and z
    initG = ones(k,1) * 0.2;
    initZ = 0.1;
    initvec = [initG; initZ];
    
    both = @(inputvec) logp_and_grad(inputvec,x,A,B,sigma_x,sigma_v,g_shape,g_scale,z_shape,z_scale);
    s = nuts_da(both,L,burnin,initvec',target_acceptance);
    gsamp = s(:,1:k);
    zsamp = s(:,end);    
end

function [logp,grad] = logp_and_grad(inputvec,x,A,B,sigma_x,sigma_v,g_shape,g_scale,z_shape,z_scale)
    inputvec = inputvec';
    logp = msmLogPostGZ(inputvec(1:end-1,1),inputvec(end,1),x,A,B,sigma_x,sigma_v,g_shape,g_scale,z_shape,z_scale);
    grad = msmLogPostGZGrad(inputvec(1:end-1,1),inputvec(end,1),x,A,B,sigma_x,sigma_v,g_shape,g_scale,z_shape,z_scale)';
end
