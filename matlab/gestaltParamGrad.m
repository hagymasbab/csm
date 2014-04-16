function grad = gestaltParamGrad(ge,samples,cholesky,varargin)
    parser = inputParser;
    addParamValue(parser,'verbose',0,@isnumeric);
    addParamValue(parser,'precision',0,@islogical);
    parse(parser,varargin{:});
    verb = parser.Results.verbose;    
    precision = parser.Results.precision;    

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
        GG = squeeze(samples(n,:,1:ge.k));
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
            for i=1:ge.k
                grad{i} = grad{i} - (g(i,1) * matr);
            end
        end
    end
    for i=1:ge.k
        grad{i} = grad{i} * cholesky{i} / L;
    end
    if verb > 0
        fprintf('\n');
    end
end
                