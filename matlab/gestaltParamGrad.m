function grad = gestaltParamGrad(ge,samples,cholesky,varargin)
    % gradient of the complete data log-likelihood with respect to the
    % cholesky decomposition of the covariance components

    parser = inputParser;
    addParamValue(parser,'verbose',0,@isnumeric);
    addParamValue(parser,'precision',false,@islogical);
    parse(parser,varargin{:});
    verb = parser.Results.verbose;    
    precision = parser.Results.precision;    

    %size(samples)
    L = size(samples,2);
    N = size(samples,1);
    grad = cell(1,ge.k);
    cc = cell(1,ge.k);
    for i=1:ge.k
        grad{i} = zeros(ge.Dv);
        cc{i} = cholesky{i}' * cholesky{i};
    end    
    
    % calculate U_hat
    U_hat = cell(ge.k,ge.Dv,ge.Dv);
    for kk=1:ge.k
        if verb > 0
            printCounter(kk,'stringVal','gradComp','maxVal',ge.k,'newLine',false);
        end
        for i=1:ge.Dv
            % this is an upper triangle matrix
            for j=i:ge.Dv
                U_hat{kk,i,j} = derivQuadByElement(cholesky{kk},i,j);
            end
        end
    end
    
    
    for n=1:N
        if verb > 0
            printCounter(n,'stringVal','gradData','maxVal',N,'newLine',false);
        end
        GG = reshape(samples(n,:,1:ge.k),L,ge.k); % squeeze doesn't work for L=1
        for l=1:L
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
                dLdC = (-1/(2*L)) * (ge.B * iCv - quad);                
            else                        
                dLdC = (1/L) * (ge.B * eye(ge.Dv)) / CvP - VV;                
            end          
            % calculate the derivative of the covariance matrix w.r.t. each
            % element of each covariance component
            for kk=1:ge.k
                for i=1:ge.Dv
                    % this is an upper triangle matrix
                    for j=i:ge.Dv
                        %U_hat = derivQuadByElement(cholesky{kk},i,j);
                        dCdu = g(kk,:) * U_hat{kk,i,j};
                        
                        % gradient of the log-gaussian w.r.t. the actual
                        % element of the actual component, summed over
                        % samples
                        grad{kk}(i,j) = grad{kk}(i,j) + trace(dLdC' * dCdu);
                    end
                end
            end
        end
    end

end
                