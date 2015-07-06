function [gsamp,zsamp] = gestaltHamiltonianGZ(x,ge,L,burnin,randseed)
    
    setrandseed(randseed);
    x = reshape(x,ge.Dv,1);    
    
    % set initial values for g and z
    initG = ones(ge.k,1) * 0.2;
    initZ = 0.1;
    initvec = [initG; initZ];
    
    both = @(inputvec) logp_and_grad(inputvec,ge,x);
    s = nuts_da(both,L,burnin,initvec');
    gsamp = s(:,1:ge.k);
    zsamp = s(:,end);    
end

function [logp,grad] = logp_and_grad(inputvec,ge,x)
    inputvec = inputvec';
    logp = gestaltLogPostGZ(inputvec(1:ge.k,1),inputvec(end,1),x,ge);
    grad = gestaltLogPostGZGrad(inputvec(1:ge.k,1),inputvec(end,1),x,ge)';
end
