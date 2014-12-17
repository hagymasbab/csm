function [vsamp,gsamp,zsamp] = gestaltHamiltonian(ge,X,nSamp,varargin)
    % TODO set X
    
    s = zeros(nSamp,ge.k + ge.B*ge.Dv + 1);
    
    % TODO set initial values for g and z
    
    % TODO set logpdf and grad functions
    
    % set bounds
    bounds = [[(1:ge.k)';ge.k+ge.B*ge.Dv+1] repmat([0 Inf],ge.k+1,1)];
    
    % TODO call sampler
    [s,rr] = hamiltonianMC(init,logpdf,grad,nSamp,lfSteps,stepSize,'bounds',bounds);
end