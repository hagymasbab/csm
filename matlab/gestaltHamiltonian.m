function [vsamp,gsamp,zsamp,rr] = gestaltHamiltonian(ge,X,nSamp,varargin)
    parser = inputParser;
    addParameter(parser,'verbose',0,@isnumeric);
    addParameter(parser,'lfSteps',20,@isnumeric);
    addParameter(parser,'stepSize',1e-4,@isnumeric);
    parse(parser,varargin{:});
    params = parser.Results;

    % set X
    if ndims(X) == 3
        X = reshape(X,ge.B,ge.Dv);
    end    
    
    % set initial values for v,g and z
    initG = ones(ge.k,1) * 0.1;
    initZ = 0.1;
    initV = ones(ge.B,ge.Dv) * 0.1;
    initvec = [initG;reshape(initV,ge.B*ge.Dv,1);initZ];
    
    % set logpdf and grad functions
    logpdf = @(inputvec) gestaltFullLogPosterior(ge,X,reshape(inputvec(ge.k+1:end-1,1),ge.B,ge.Dv),inputvec(1:ge.k,1),inputvec(end,1),[]);
    grad = @(inputvec) gestaltFullLogPosteriorGrad(ge,X,reshape(inputvec(ge.k+1:end-1,1),ge.B,ge.Dv),inputvec(1:ge.k,1),inputvec(end,1),[]);
    
    % set bounds
    bounds = [[(1:ge.k)';ge.k+ge.B*ge.Dv+1] repmat([0 Inf],ge.k+1,1)];
    
    % call sampler
    [s,rr] = hamiltonianMC(initvec,logpdf,grad,nSamp,params.lfSteps,params.stepSize,'bounds',bounds,'verbose',params.verbose);
    gsamp = s(:,1:ge.k);
    vsamp = reshape(s(:,ge.k+1:end-1),[nSamp ge.B ge.Dv]);
    zsamp = s(:,end);
end
