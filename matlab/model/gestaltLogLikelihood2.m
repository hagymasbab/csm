function ll = gestaltLogLikelihood2(ge,L,data,cholesky,varargin)    
    parser = inputParser;   
    addParameter(parser,'verbose',0,@isnumeric);  
    addParameter(parser,'loadSamples',false,@islogical);
    addParameter(parser,'randseed','leave'); 
    addParameter(parser,'scientific','true'); % dummy parameter
    addParameter(parser,'method','intuition'); % TODO this stuff
    parse(parser,varargin{:});        
    params = parser.Results;  
    
    setrandseed(params.randseed);
    N = size(data,1);
    if ge.B ~= 1 || ge.nullComponent
        error('not implemented');
    end
    data = reshape(data,N,ge.Dx);
    
    cc = cell(1,ge.k);    
    for i=1:ge.k
        cc{i} = cholesky{i}' * cholesky{i};
    end  
    
    logConstant = -N * log(L);
    
    if params.loadSamples
        load('prior_samples.mat');
    else
        G = gestaltSamplePriorG(ge,L);
        Z = gamrnd(ge.z_shape,ge.z_scale,[L 1]);
        save('bin/prior_samples.mat','G','Z');
    end
    
    inverse_covariances = cell(1,L);
    plogdet = zeros(L,1);
    loghz = ones(L,1);
    for l=1:L
        cv = componentSum(G(l,:)',cc);    
        c_right = Z(l,1)^2 * (ge.A * cv * ge.A');
        c_right = (c_right + c_right') / 2;
        nCl = ge.obsVar * eye(ge.Dv) + c_right;            
        Cl = nearestSPD(nCl);
        plogdet(l) = stableLogdet(Cl,'scaling','unknown');
        inverse_covariances{l} = stableInverse(Cl);
    end
    
    ll = 0;
    for n = 1:N
        if params.verbose == 1
            printCounter(n,'maxVal',N,'stringVal','Likelihood datapoint');
        end
        
        ll_part_coeff = 0;
        ll_part_expo = 0;
        x_n = data(n,:)';
        
        for l = 1:L                                  
            [ll_act_coeff,ll_act_expo] = stableMvnpdf(x_n,zeros(ge.Dx,1),inverse_covariances{l},'invertedC',true,'scientific',true,'precompLogdet',plogdet(l));
            [ll_part_coeff,ll_part_expo] = sumSciNot(ll_part_coeff,ll_part_expo,ll_act_coeff,ll_act_expo);
        end
        ll = ll + scinot2log(ll_part_coeff,ll_part_expo);
    end
    
    ll = ll + logConstant;
end