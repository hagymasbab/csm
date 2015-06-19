function grad = gestaltLogLikelihoodGradient(ge,L,data,cholesky,varargin)
    parser = inputParser;   
    addParameter(parser,'verbose',0,@isnumeric);  
    addParameter(parser,'loadSamples',false,@islogical);
    addParameter(parser,'randseed','leave');      
    addParameter(parser,'method','scinot');      % dummy parameter
    addParameter(parser,'template',true(ge.Dv));
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
    
    pA = pinv(ge.A);
    ATA = ge.A' * ge.A;     
    siATA = ge.obsVar * stableInverse(ATA);          
    
    if params.loadSamples
        load('prior_samples.mat');
    else
        G = gestaltSamplePriorG(ge,L,'checkValues',false);
        Z = gamrnd(ge.z_shape,ge.z_scale,[L 1]);
        save('bin/prior_samples.mat','G','Z');
    end
               
    % stuff that only depends on the samples, not the data
    loghz = zeros(L,1);
    plogdet = zeros(L,1);
    inverse_covariances = cell(1,L);
    for l=1:L              
        if params.verbose == 1
            printCounter(l,'stringVal','Sample','maxVal',L,'newLine',false);
        end
        g_l = G(l,:)';
        z_l = Z(l,1);
        cv = componentSum(g_l,cc);   
        leftmat = siATA / Z(l)^2;
        nCl = leftmat + cv;                    
        Cl = nearestSPD(nCl);        
        plogdet(l) = stableLogdet(Cl,'scaling','up');
        inverse_covariances{l} = stableInverse(Cl);
        loghz(l) = -ge.Dv * log(z_l);
    end    
    
    
    % data-dependent stuff
    M_cell = cell(1,ge.k);
    for kk = 1:ge.k
        M_cell{kk} = zeros(ge.Dv);
    end
    h_times_N_coeff = zeros(N,L);        
    h_times_N_expo = zeros(N,L);
    for n = 1:N
        if params.verbose == 1
            printCounter(n,'stringVal','Datapoint','maxVal',N,'newLine',false);
        end
        x_n = data(n,:)';
        pAx = pA * x_n;
        
        Li_coeff = 0;
        Li_expo = 0;
        for l = 1:L
            f_l = pAx / Z(l);
            iC_l = inverse_covariances{l};
            [coeff_n,expo_n] = stableMvnpdf(f_l,zeros(ge.Dv,1),iC_l,'scientific',true,'invertedC',true,'precompLogdet',plogdet(l));
            [coeff_h,expo_h] = sciNot(loghz(l),true);   
            [h_times_N_coeff(n,l),h_times_N_expo(n,l)] = prodSciNot([coeff_n coeff_h],[expo_n expo_h]);
            [Li_coeff,Li_expo] = sumSciNot(Li_coeff,Li_expo,h_times_N_coeff(n,l),h_times_N_expo(n,l));
        end        
        
        M_part_cell = cell(1,ge.k);
        for kkk = 1:ge.k
            M_part_cell{kkk} = zeros(ge.Dv);
        end
                
        for l = 1:L      
            
            f_l = pAx / Z(l);
            iC_l = inverse_covariances{l};
            iCf = iC_l * f_l;       
            
            [hNperL_coeff,hNperL_expo] = prodSciNot([h_times_N_coeff(n,l) 1/Li_coeff],[h_times_N_expo(n,l) -Li_expo]);
            scalar_term = hNperL_coeff * 10^hNperL_expo;
                       
            matrix_term = iC_l - iCf * iCf';
            s_times_M = scalar_term * matrix_term;
            
            for kkk = 1:ge.k
                M_part_update = G(l,kkk) * s_times_M;
                M_part_cell{kkk} = M_part_cell{kkk} + M_part_update;
            end                                                                
        end
        
        M_cell = celladd(M_cell,1,M_part_cell,-1/2);
    end
        
    % cycle over all the individual parameters
    grad = cell(1,ge.k);          
    grad{kk} = zeros(ge.Dv);
    for kk = 1:ge.k     
        grad{kk} = zeros(ge.Dv);
    end
    for i = 1:ge.Dv
        if params.verbose == 1
            printCounter(i,'stringVal','Row','maxVal',ge.Dv,'newLine',true);
        end
        for j = i:ge.Dv      
            if params.template(i,j)                               
                for kk = 1:ge.k     
                    %M_cell{kk}(j,:)
                    grad{kk}(i,j) = sum(cholesky{kk}(i,:) .* M_cell{kk}(j,:) + cholesky{kk}(i,:) .* M_cell{kk}(:,j)');
                end
            end
        end
    end    
end