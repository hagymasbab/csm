function grad = gestaltParamGrad(ge,vsamp,gsamp,cholesky,varargin)
    % gradient of the complete data log-likelihood with respect to the
    % cholesky decomposition of the covariance components

    parser = inputParser;
    addParamValue(parser,'verbose',0,@isnumeric);
    addParamValue(parser,'precision',false,@islogical);
    parse(parser,varargin{:});
    verb = parser.Results.verbose;    
    precision = parser.Results.precision;    
                
    L = size(vsamp,2);
    N = size(vsamp,1);
        
    cc = cell(1,ge.k);    
    if ge.nullComponent
        effective_k = ge.k-1;
        cc{end} = ge.cc{end};
    else
        effective_k = ge.k;
    end
    
    grad = cell(1,effective_k);
    dLdC = cell(1,effective_k);
    for i=1:effective_k
        grad{i} = zeros(ge.Dv);
        cc{i} = cholesky{i}' * cholesky{i};
        dLdC{i} = zeros(ge.Dv);
    end            
    
    for n=1:N        

        for l=1:L
            if verb > 0
                printCounter((n-1)*L+l,'stringVal','gradSample','maxVal',N*L,'newLine',false);
            end

            g = reshape(gsamp(n,l,:),ge.k,1);
            V = reshape(vsamp(n,l,:,:),ge.B,ge.Dv);
            
            % calculate \sum_{b=1}^B v^{m,b} v^{(m,b)T}            
            VV = V'*V;
            
            % calculate the covariance matrix
            CvP = componentSum(g,cc);
            
            if ~precision                
                % derivative of the log-gaussian formula w.r.t the
                % covariance matrix
                % dLdC = (-1/(2*L)) * ( (ge.B * eye(ge.Dv)) / CvP - (CvP \ VV) / CvP );
                
                iCv = inv(CvP);
                quad = 0;
                for b=1:ge.B
                    v = V(b,:)';
                    quad = quad + iCv * (v * v') * iCv;
                end                
                actmat = (ge.B * iCv - quad);
                for kk = 1:effective_k
                    dLdC{kk} = dLdC{kk} + g(kk,:) * actmat;               
                end
            else    
                % TODO replace this with writing to stored matrices
                dLdC = dLdC + (1/L) * (ge.B * eye(ge.Dv)) / CvP - VV;                
            end      
        end
    end
        
    % calculate the derivative of the covariance matrix w.r.t. each
    % element of each covariance component
    parfor kk=1:effective_k
        if verb > 0
            printCounter(kk,'stringVal','gradComp','maxVal',effective_k,'newLine',true);
        end
        for i=1:ge.Dv
            % this is an upper triangle matrix
            for j=i:ge.Dv
                % the computation above is equivalent with this much faster one
                grad{kk}(i,j) = sum(cholesky{kk}(i,:) .* dLdC{kk}(j,:) + cholesky{kk}(i,:) .* dLdC{kk}(:,j)');
            end
        end
        grad{kk} = (-1/(2*L)) * grad{kk};
    end
            
end
                