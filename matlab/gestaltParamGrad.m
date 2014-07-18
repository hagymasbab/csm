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
    if verb > 0
        fprintf('Sample %d/', N*L);
    end
    for n=1:N
        if verb > 0
            printCounter(n);
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
                dLdC = (-1/(2*L)) * ( (ge.B * eye(ge.Dv)) / CvP - (CvP \ VV) / CvP );
                
                %iCv = inv(CvP);
                %matr = ge.B * iCv - iCv * VV * iCv;                
            else                        
                dLdC = (ge.B * eye(ge.Dv)) / CvP - VV;                
            end          
            % calculate the derivative of the covariance matrix w.r.t. each
            % element of each covariance component
            for kk=i:ge.k
                for i=1:ge.Dv
                    for j=1:ge.Dv
                        U_hat = derivQuadByElement(cholesky{kk},i,j);
                        dCdu = grad{kk}(i,j) * U_hat;
                        
                        % gradient of the log-gaussian w.r.t. the actual
                        % element of the actual component, summed over
                        % samples
                        grad{kk}(i,j) = grad{kk}(i,j) + trace(dLdC' * dCdu);
                    end
                end
            end
        end
    end

    if verb > 0
        fprintf('\n');
    end
end
                