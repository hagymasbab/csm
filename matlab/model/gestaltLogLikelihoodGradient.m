function grad = gestaltLogLikelihoodGradient(ge,L,data,cholesky,varargin)
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
        G = gestaltSamplePriorG(ge,L);
        Z = gamrnd(ge.z_shape,ge.z_scale,[L 1]);
        save('bin/prior_samples.mat','G','Z');
    end
    Z2 = Z .^ 2;       
    
    covariances = zeros(L,ge.Dv,ge.Dv);
    inverse_covariances = zeros(L,ge.Dv,ge.Dv);
    hz = zeros(L,1);
    parfor l=1:L
        g_l = G(l,:)';
        z_l = Z(l,1);
        cv = componentSum(g_l,cc);        
        leftmat = siATA / Z2(l);
        nCl = leftmat + cv;                    
        Cl = nearestSPD(nCl);        
        [~,e1] = cholcov(leftmat);
        [~,e2] = cholcov(cv);
        [~,e3] = cholcov(nCl);
        [~,e4] = cholcov(Cl);
        fprintf('Condition siATA/z %e Cv %e nCl %e Cl %e\n',rcond(leftmat),rcond(cv),rcond(nCl),rcond(Cl));
        fprintf('Posdef siATA/z %d Cv %d nCl %d Cl %d\n',e1,e2,e3,e4);
        covariances(l,:,:) = Cl;
        inverse_covariances(l,:,:) = stableInverse(reshape(covariances(l,:,:),ge.Dv,ge.Dv));
        %inverse_covariances(l,:,:) = inv(reshape(covariances(l,:,:),ge.Dv,ge.Dv));
        %hz(l) = idATA / (z_l^ge.Dv);
        hz(l) = 1 / (z_l^ge.Dv);
    end    
    if any(isnan(hz))
        error('NaN in hz');
    end
    %Z.^ge.Dv
    %Z
    %hz
        
    M = zeros(ge.k,ge.Dv,ge.Dv);
    parfor n = 1:N
        if params.verbose == 1
            printCounter(n,'stringVal','Datapoint','maxVal',N);
        end
        x_n = data(n,:)';
        pAx = pA * x_n;
        Li_n = 0;
        M_part = zeros(ge.k,ge.Dv,ge.Dv);
        zero_scalars = 0;
        for l = 1:L      
            %printCounter((n-1)*L+l,'maxVal',N*L,'stringVal','e')
            f_l = pAx / Z(l);
            C_l = reshape(covariances(l,:,:),ge.Dv,ge.Dv);
            iC_l = reshape(inverse_covariances(l,:,:),ge.Dv,ge.Dv);
            iCf = iC_l * f_l;
            
%             size(C_l)
%             issymmetric(C_l)
%              [~,e] = cholcov(C_l);
%              e
%              [~,e] = cholcov(iC_l);
%              e                                            
            
            N_f = mvnpdf(f_l',zeros(1,ge.Dv),C_l);
            [a,e] = stableMvnpdf(f_l,zeros(ge.Dv,1),iC_l,true,true);            
            %e
            
            if any(isnan(N_f))
                error('NaN in N_f,n %d %d',n,l);
            end
            %N_f = stableMvnpdf(f_l,zeros(ge.Dv,1),C_l,false,false)
            if N_f == 0                
                scalar_term = 0;
                zero_scalars = zero_scalars + 1;
            else                
                scalar_term = hz(l) * N_f;   
            end
            if any(isnan(scalar_term))
                hz(l)
                N_f
                error('NaN in scalar,n %d %d',n,l);
            end
            matrix_term = iC_l - iCf * iCf';
            
            for kk = 1:ge.k
                M_part(kk,:,:) = M_part(kk,:,:) + reshape(G(l,kk) * scalar_term * matrix_term,1,ge.Dv,ge.Dv);
            end
            
            Li_n = Li_n + scalar_term;            
        end
        zero_scalars
        % this is as approximation
        if Li_n ~= 0
            M = M + M_part / Li_n;      
        end
    end
    
    if any(isnan(M))
        error('NaN in M');
    end
    
    M = -M / 2;
    
%     viewImageSet(M)
%     pause
    
    grad = cell(1,ge.k);
    for kk = 1:ge.k
        grad{kk} = zeros(ge.Dv);
        M_k = reshape(M(kk,:,:),ge.Dv,ge.Dv);
        for i = 1:ge.Dv
            for j = i:ge.Dv                
                %U_hat = derivQuadByElement(cholesky{k},i,j);
                %grad{kk}(i,j) = trace(M*U_hat);
                % this should be equivalent with Tr(M U_hat)                
                grad{kk}(i,j) = sum(cholesky{kk}(i,:) .* M_k(j,:) + cholesky{kk}(i,:) .* M_k(:,j)');
            end
        end
    end
    
    if any(isnan(cell2mat(grad)))
        error('NaN in grad');
    end
end