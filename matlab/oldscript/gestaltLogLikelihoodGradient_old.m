function grad = gestaltLogLikelihoodGradient(ge,L,data,cholesky,varargin)
    parser = inputParser;   
    addParameter(parser,'verbose',0,@isnumeric);  
    addParameter(parser,'loadSamples',false,@islogical);
    addParameter(parser,'randseed','leave');      
    addParameter(parser,'method','scinot');      
    addParameter(parser,'template',[]);
    parse(parser,varargin{:});        
    params = parser.Results;  
    
    setrandseed(params.randseed);
    N = size(data,1);
    if ge.B ~= 1 || ge.nullComponent
        error('not implemented');
    end
    data = reshape(data,N,ge.Dx);
    
    if isempty(params.template)
        params.template = true(ge.Dv);
    end
    paramNum = ge.k + sum(params.template(:));
    
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
    Z2 = Z .^ 2;       
    
    covariances = zeros(L,ge.Dv,ge.Dv);
    inverse_covariances = zeros(L,ge.Dv,ge.Dv);
    % TODO csak loghz kell
    hz = zeros(L,1);
    loghz = zeros(L,1);
    for l=1:L              
        if params.verbose == 1
            printCounter(l,'stringVal','Sample','maxVal',L,'newLine',false);
        end
        g_l = G(l,:)';
        z_l = Z(l,1);
        cv = componentSum(g_l,cc);        
        if sum(sum(abs(cv)>0)) < 0.2 * length(cv);
            cv = sparse(cv);
        end
        leftmat = siATA / Z2(l);
        nCl = leftmat + cv;                    
        Cl = nearestSPD(nCl);        
%         [~,e1] = cholcov(leftmat);
%         [~,e2] = cholcov(cv);
%         [~,e3] = cholcov(nCl);
%         [~,e4] = cholcov(Cl);
%         fprintf('Condition siATA/z %e Cv %e nCl %e Cl %e\n',rcond(leftmat),rcond(cv),rcond(nCl),rcond(Cl));
%         fprintf('Posdef siATA/z %d Cv %d nCl %d Cl %d\n',e1,e2,e3,e4);
        covariances(l,:,:) = Cl;
        inverse_covariances(l,:,:) = stableInverse(reshape(covariances(l,:,:),ge.Dv,ge.Dv));
        %inverse_covariances(l,:,:) = inv(reshape(covariances(l,:,:),ge.Dv,ge.Dv));
        %hz(l) = idATA / (z_l^ge.Dv);
        hz(l) = 1 / (z_l^ge.Dv);
        loghz(l) = -ge.Dv * log(z_l);
    end    
    if any(isnan(hz))
        error('NaN in hz');
    end
    %Z.^ge.Dv
    %Z
    %loghz
        
    M = zeros(ge.k,ge.Dv,ge.Dv);
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
        
        if strcmp(params.method,'scinot')
            % first calculate h_l * N_l for numerical stability
%             h_times_N_coeff = zeros(L,1);        
%             h_times_N_expo = zeros(L,1);        
            Li_coeff = 0;
            Li_expo = 0;
            for l = 1:L
                f_l = pAx / Z(l);
                iC_l = reshape(inverse_covariances(l,:,:),ge.Dv,ge.Dv);
                [coeff_n,expo_n] = stableMvnpdf(f_l,zeros(ge.Dv,1),iC_l,true,true);  
                %expo_n
                %[coeff_h,expo_h] = sciNot(hz(l),false);
                [coeff_h,expo_h] = sciNot(loghz(l),true);   
                %expo_h
                [h_times_N_coeff(n,l),h_times_N_expo(n,l)] = prodSciNot([coeff_n coeff_h],[expo_n expo_h]);
                [Li_coeff,Li_expo] = sumSciNot(Li_coeff,Li_expo,h_times_N_coeff(n,l),h_times_N_expo(n,l));
            end        
    %         h_times_N_coeff
    %         h_times_N_expo
    %         
    %         Li_coeff
        else
            Li_n = 0;        
        end
        
%         M_part = zeros(ge.k,ge.Dv,ge.Dv);
        M_part_cell = cell(1,ge.k);
        for kk = 1:ge.k
            M_part_cell{kk} = zeros(ge.Dv);
        end
        zero_scalars = 0;
                
        for l = 1:L      
            
            f_l = pAx / Z(l);
            iC_l = reshape(inverse_covariances(l,:,:),ge.Dv,ge.Dv);
            iCf = iC_l * f_l;       
            
            if strcmp(params.method,'scinot')
                [hNperL_coeff,hNperL_expo] = prodSciNot([h_times_N_coeff(n,l) 1/Li_coeff],[h_times_N_expo(n,l) -Li_expo]);
                scalar_term = hNperL_coeff * 10^hNperL_expo;
            else
                C_l = reshape(covariances(l,:,:),ge.Dv,ge.Dv);
                N_f = mvnpdf(f_l',zeros(1,ge.Dv),C_l);                

                if any(isnan(N_f))
                    error('NaN in N_f,n %d %d',n,l);
                end
                
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
                
                 Li_n = Li_n + scalar_term;  
            end
            
            matrix_term = iC_l - iCf * iCf';
            
            for kk = 1:ge.k
                M_part_update = G(l,kk) * scalar_term * matrix_term;
%                 M_part_update2 = reshape(M_part_update,1,ge.Dv,ge.Dv);
%                 M_part(kk,:,:) = M_part(kk,:,:) + M_part_update2;
                M_part_cell{kk} = M_part_cell{kk} + M_part_update;
            end                                                                
        end
        
        %zero_scalars
        % this is as approximation
        if ~strcmp(params.method,'scinot') && Li_n ~= 0
%             M_part = M_part / Li_n;      
            celladd(M_cell,1,M_part_cell,-1/(2*Li_n));
        else        
%             M = M + M_part;
            celladd(M_cell,1,M_part_cell,-1/2);
        end
    end
    
    if any(isnan(cell2mat(M_cell)))
        error('NaN in M');
    end
    
%     M = -M / 2;
    
%     viewImageSet(M)
%     pause
    
%     counter = 1;
%     grad = cell(1,ge.k);
%     for kk = 1:ge.k        
%         if params.verbose == 1
%             printCounter(kk,'stringVal','Component','maxVal',ge.k,'newLine',true);
%         end
%         grad{kk} = zeros(ge.Dv);
% %         M_k = reshape(M(kk,:,:),ge.Dv,ge.Dv);
%         M_k = M_cell{kk};
%         for i = 1:ge.Dv
%             for j = i:ge.Dv      
%                 if params.template(i,j)                    
%                     counter = counter + 1;
%                     %U_hat = derivQuadByElement(cholesky{k},i,j);
%                     %grad{kk}(i,j) = trace(M*U_hat);
%                     % this should be equivalent with Tr(M U_hat)                
%                     grad{kk}(i,j) = sum(cholesky{kk}(i,:) .* M_k(j,:) + cholesky{kk}(i,:) .* M_k(:,j)');
%                 end
%             end
%         end
%     end
    
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
                    grad{kk}(i,j) = sum(cholesky{kk}(i,:) .* M_cell{kk}(j,:) + cholesky{kk}(i,:) .* M_cell{kk}(:,j)');
                end
            end
        end
    end
        
    if any(isnan(cell2mat(grad)))
        error('NaN in grad');
    end
end