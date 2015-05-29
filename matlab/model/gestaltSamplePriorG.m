function G = gestaltSamplePriorG(ge,L,varargin)
    parser = inputParser;
    addParameter(parser,'sampleRetry',10,@isnumeric);
    addParameter(parser,'precision',false,@islogical);
    addParameter(parser,'checkValues',true,@islogical);
    parse(parser,varargin{:});
    params = parser.Results;
    
    G = zeros(L,ge.k);
    for i=1:L
        valid = false;
        tries = 0;
        while ~valid
            if strcmp(ge.prior,'gamma')
                g_act = zeros(ge.k,1);
                for d = 1:ge.k            
                    if ge.nullComponent && d == ge.k
                        g_act(d,1) = gamrnd(ge.null_shape,ge.null_scale);
                    else
                        g_act(d,1) = gamrnd(ge.g_shape,ge.g_scale);
                    end
                end
            elseif strcmp(ge.prior,'dirichlet')
                 g_act = symmetricDirichlet(ge.sparsity,ge.k,1)';
            else
                error('Prior distribution %s not implemented',ge.prior);
            end

            if params.checkValues
                valid = gestaltCheckG(g_act,ge,params.precision);
            else
                valid = true;
            end
            tries = tries + 1;
            % if we cannot find a valid g, return an error code
            if tries > params.sampleRetry
                throw(MException('Gestalt:Prior:TooManyTries',sprintf('Number of tries to draw a valid g vector from the prior exceeded %d',params.sampleRetry)));            
            end
        end
        G(i,:) = g_act';
    end
end