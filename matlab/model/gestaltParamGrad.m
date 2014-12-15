function grad = gestaltParamGrad(ge,samples,cholesky,varargin)
    % gradient of the complete data log-likelihood with respect to the
    % cholesky decomposition of the covariance components

    parser = inputParser;
    addParamValue(parser,'verbose',0,@isnumeric);
    addParamValue(parser,'precision',false,@islogical);
    parse(parser,varargin{:});
    verb = parser.Results.verbose;    
    precision = parser.Results.precision;    
        
    if ge.nullComponent
        effective_k = ge.k-1;
    else
        effective_k = ge.k;
    end
    
    %size(samples)
    L = size(samples,2);
    N = size(samples,1);
    grad = cell(1,effective_k);
    cc = cell(1,ge.k);
    if ge.nullComponent
        cc{end} = ge.cc{end};
    end
    for i=1:effective_k
        grad{i} = zeros(ge.Dv);
        cc{i} = cholesky{i}' * cholesky{i};
    end    
    
    % calculate U_hat
%     U_hat = cell(ge.k,ge.Dv,ge.Dv);    
%     for kk=1:ge.k
%         if verb > 0
%             printCounter(kk,'stringVal','gradComp','maxVal',ge.k,'newLine',false);
%         end
%         whos
%         for i=1:ge.Dv
%             % this is an upper triangle matrix
%             for j=i:ge.Dv
%                 U_hat{kk,i,j} = derivQuadByElement(cholesky{kk},i,j);
%             end
%         end
%     end
    
    dLdC = cell(1,effective_k);
    for kk = 1:effective_k
        dLdC{kk} = zeros(ge.Dv);
    end
    
    for n=1:N        
        GG = reshape(samples(n,:,1:ge.k),L,ge.k); % squeeze doesn't work for L=1
        parfor l=1:L
            if verb > 0
                printCounter((n-1)*L+l,'stringVal','gradSample','maxVal',N*L,'newLine',false);
            end
            % retrieve g^m
            g = GG(l,:)';
            % calculate \sum_{b=1}^B v^{m,b} v^{(m,b)T}
            V = reshape(samples(n,l,ge.k+1:ge.k+ge.Dv*ge.B),ge.B,ge.Dv);
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
                    %size(v)
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
%                 U_hat = derivQuadByElement(cholesky{kk},i,j);
%                 %dCdu = U_hat{kk,i,j};
% 
%                 % gradient of the log-gaussian w.r.t. the actual
%                 % element of the actual component, summed over
%                 % samples
%                 grad{kk}(i,j) = trace(dLdC{kk}' * U_hat);

                % the computation above is equivalent with this much faster one
                grad{kk}(i,j) = sum(cholesky{kk}(i,:) .* dLdC{kk}(j,:) + cholesky{kk}(i,:) .* dLdC{kk}(:,j)');
            end
        end
        grad{kk} = (-1/(2*L)) * grad{kk};
    end
            
end
                