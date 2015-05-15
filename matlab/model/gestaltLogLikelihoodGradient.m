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
    
    if params.loadSamples
        load('prior_samples.mat');
    else
        G = gestaltSamplePriorG(ge,L);
        Z = gamrnd(ge.z_shape,ge.z_scale,[L 1]);
        save('bin/prior_samples.mat','G','Z');
    end
    Z_2 = Z .^ 2;
    
    pA = pinv(ge.A);
    ATA = ge.A' * ge.A;
    siATA = ge.obsVar * inv(ATA);
    idATA = 1 / sqrt( det(ATA) );
    
    % TODO precompute h_l, C_l and iC_l    
    
    M = zeros(ge.k,ge.Dv,ge.Dv);
    for n = 1:N
        if params.verbose == 1
            printCounter(n,'stringVal','Datapoint','maxVal',N);
        end
        x_n = data(n,:)';
        pAx = pA * x_n;
        Li_n = 0;
        M_part = zeros(ge.k,ge.Dv,ge.Dv);
        for l = 1:L            
            g_l = G(l,:)';
            z_l = Z(l,1);
            z_l2 = Z_2(l,1);
            
            h_l = 1 / (z_l^ge.Dv * idATA);
            f_l = pAx / z_l;
            C_v = componentSum(g_l,cc);
            C_l = siATA / z_l2 + C_v;
            
%             [~,e] = chol(C_l);
%             if e ~= 0
%                 C_l = nearestSPD(C_l);
%             end
                        
            iC_l = inv(C_l);
            iCf = iC_l * f_l;
            
            N_f = mvnpdf(f_l',zeros(1,ge.Dv),C_l);
            scalar_term = h_l * N_f;            
            
            matrix_term = iC_l - iCf * iCf';
            
            %size(matrix_term)
            %size(M_part(1,:,:))
            
            for kk = 1:ge.k
                M_part(kk,:,:) = M_part(kk,:,:) + reshape(g_l(kk,1) * scalar_term * matrix_term,1,ge.Dv,ge.Dv);
            end
            
            Li_n = Li_n + scalar_term;            
        end
        
        M = M + M_part / Li_n;        
    end
    
    M = -M / 2;
    
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
end