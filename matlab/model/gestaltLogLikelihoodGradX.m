function grad = gestaltLogLikelihoodGradX(ge,L,data,cholesky,sigma_x,varargin)
    parser = inputParser;   
    addParameter(parser,'verbose',0,@isnumeric);  
    addParameter(parser,'loadSamples',false,@islogical);
    addParameter(parser,'randseed','leave');              
    parse(parser,varargin{:});        
    params = parser.Results;  
    
    setrandseed(params.randseed);
    N = size(data,1);
    if ge.B ~= 1 || ge.nullComponent
        error('not implemented');
    end
    data = reshape(data,N,ge.Dx);
    
    cc = cholcell(cholesky);                       
    
    if params.loadSamples
        load('prior_samples.mat');
    else
        G = gestaltSamplePriorG(ge,L,'checkValues',false);
        Z = gamrnd(ge.z_shape,ge.z_scale,[L 1]);
        save('bin/prior_samples.mat','G','Z');
    end
               
    % stuff that only depends on the samples, not the data    
    plogdet = zeros(L,1);
    inverse_covariances = cell(1,L);
    for l=1:L              
        if params.verbose == 1
            printCounter(l,'stringVal','Sample','maxVal',L,'newLine',false);
        end
        g_l = G(l,:)';
        z_l = Z(l,1);
        cv = componentSum(g_l,cc);   
        leftmat = sigma_x * eye(ge.Dv);
        nCl = leftmat + z_l^2 * (ge.A * cv * ge.A');                    
        Cl = nearestSPD(nCl);        
        plogdet(l) = stableLogdet(Cl,'scaling','up');
        inverse_covariances{l} = stableInverse(Cl);        
    end    
        
    grad = 0;    
    % data-dependent stuff
    for n = 1:N
        if params.verbose == 1
            printCounter(n,'stringVal','Datapoint','maxVal',N,'newLine',false);
        end
        x_n = data(n,:)';
        xxT = x_n * x_n';
        
        num_coeff = 0;
        num_expo = 0;
        den_coeff = 0;
        den_expo = 0;

        for l = 1:L
            iC_l = inverse_covariances{l};
            [coeff_n,expo_n] = stableMvnpdf(x_n,zeros(ge.Dv,1),iC_l,'scientific',true,'invertedC',true,'precompLogdet',plogdet(l));           
            [den_coeff,den_expo] = sumSciNot(den_coeff,den_expo,coeff_n,expo_n);
            
            matrix_trace = trace(iC_l + iC_l*xxT*iC_l);
            [trace_coeff,trace_expo] = sciNot(matrix_trace,false);
            [ntr_coeff,ntr_expo] = prodSciNot([coeff_n trace_coeff],[expo_n trace_expo]);
            [num_coeff,num_expo] = sumSciNot(num_coeff,num_expo,ntr_coeff,ntr_expo);            
        end        
        
        [point_coeff,point_expo] = prodSciNot([num_coeff,1/den_coeff],[num_expo,-den_expo]);
        grad = grad + point_coeff*10^point_expo;
    end
        
    grad = -grad/2;
  
end