function ll = gestaltLogLikelihood2(ge,L,data,cholesky,varargin)    
    parser = inputParser;   
    addParameter(parser,'verbose',0,@isnumeric);  
    addParameter(parser,'loadSamples',false,@islogical);
    addParameter(parser,'randseed','leave'); 
    addParameter(parser,'scientific','true'); % dummy parameter
    addParameter(parser,'method','intuition'); 
    addParameter(parser,'sigma',0,@isnumeric); 
    parse(parser,varargin{:});        
    params = parser.Results;  
    
    setrandseed(params.randseed);
    N = size(data,1);
    if ge.B ~= 1 || ge.nullComponent
        error('not implemented');
    end
    data = reshape(data,N,ge.Dx);
    
    if params.sigma > 0
        ge.obsVar = params.sigma;
    end
    
    cc = cholcell(cholesky); 
    
    pA = pinv(ge.A);
    ATA = ge.A' * ge.A;     
    siATA = ge.obsVar * stableInverse(ATA);          
    
    if strcmp(params.method,'intuition')
        logConstant = -N * log(L);
    else
        logConstant = -N * (L + stableLogdet(ATA)/2);
    end
    
    if params.loadSamples
        load('prior_samples.mat');
    else
        G = gestaltSamplePriorG(ge,L,'checkValues',false);
        Z = gamrnd(ge.z_shape,ge.z_scale,[L 1]);
        save('bin/prior_samples.mat','G','Z');
    end
    
    inverse_covariances = cell(1,L);
    plogdet = zeros(L,1);
    loghz = zeros(L,1);
    for l=1:L
        cv = componentSum(G(l,:)',cc);    
        if strcmp(params.method,'intuition')
            c_right = Z(l,1)^2 * (ge.A * cv * ge.A');
            c_right = (c_right + c_right') / 2;
            nCl = ge.obsVar * eye(ge.Dv) + c_right;            
        else
            c_left = siATA / Z(l,1)^2;
            nCl = c_left + cv;
            loghz(l) = -ge.Dv * log(Z(l,1));
        end       
        Cl = nearestSPD(nCl);
%         log10(rcond(Cl))
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
        
        if ~strcmp(params.method,'intuition')
            x_n = pA * x_n;
        end
        
        for l = 1:L           
            if ~strcmp(params.method,'intuition')
                x_n = x_n / Z(l,1);
            end
            [coeff_n,expo_n] = stableMvnpdf(x_n,zeros(ge.Dx,1),inverse_covariances{l},'invertedC',true,'scientific',true,'precompLogdet',plogdet(l));
            [coeff_h,expo_h] = sciNot(loghz(l),true);   
            [ll_act_coeff,ll_act_expo] = prodSciNot([coeff_n coeff_h],[expo_n expo_h]);            
            [ll_part_coeff,ll_part_expo] = sumSciNot(ll_part_coeff,ll_part_expo,ll_act_coeff,ll_act_expo);
        end
        ll = ll + scinot2log(ll_part_coeff,ll_part_expo);
    end
    
    ll = ll + logConstant;
end