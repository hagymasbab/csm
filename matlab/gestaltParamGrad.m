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
            g = GG(l,:)';            
            V = reshape(samples(n,l,ge.k+1:ge.k+ge.Dv*ge.B),ge.B,ge.Dv);
            VV = V'*V;
            CvP = componentSum(g,cc);
            if ~precision                
                matr = (ge.B * eye(ge.Dv)) / CvP - (CvP \ VV) / CvP;
                %iCv = inv(CvP);
                %matr = ge.B * iCv - iCv * VV * iCv;                
            else                        
                matr = (ge.B * eye(ge.Dv)) / CvP - VV;
                
            end            
            trmatr = trace(matr);
            for kk=i:ge.k
                %grad{kk} = grad{kk} - (g(kk,1) * matr);
                for i=1:ge.Dv
                    for j=1:ge.Dv
                        grad{kk}(i,j) = grad{kk}(i,j) + g(kk,1) * trmatr;
                    end
                end
            end
        end
    end
    for kk=1:ge.k
        grad{kk} = - grad{kk} .* cholesky{kk} / L;
        %grad{kk} = grad{kk} .* cholesky{kk} / L;        
    end
    if verb > 0
        fprintf('\n');
    end
end
                