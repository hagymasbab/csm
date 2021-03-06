function ll = gestaltCompleteDataLogLikelihood(ge,samples,cholesky,varargin)
    
    parser = inputParser;
    addParamValue(parser,'verbose',0,@isnumeric);
    addParamValue(parser,'precision',false,@islogical);
    parse(parser,varargin{:});
    verb = parser.Results.verbose;    
    precision = parser.Results.precision;    

    L = size(samples,2);
    N = size(samples,1);
    ll = 0;
    cc = cell(1,ge.k);
    for i=1:ge.k
        cc{i} = cholesky{i}' * cholesky{i};
    end
    
    for n=1:N
        GG = reshape((samples(n,:,1:ge.k)),L,ge.k);
        for l=1:L
            g = GG(l,:)';
            V = reshape(samples(n,l,ge.k+1:ge.k+ge.Dv*ge.B),ge.B,ge.Dv);
            Cv = componentSum(g,cc);
            iCv = inv(Cv);
            
            quad = 0;
            for b=1:ge.B
                vb = V(b,:)';
                %size(vb)
                quad = quad + vb' * iCv * vb;
            end
            ll = ll + ge.B * log(det(Cv)) + quad;
        end
    end
    ll = (-1 / (2*L)) * ll;
end