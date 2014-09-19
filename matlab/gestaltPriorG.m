function g_act = gestaltPriorG(ge,distribution,varargin)
    parser = inputParser;
    addParamValue(parser,'sampleRetry',0,@isnumeric);
    addParamValue(parser,'precision',false,@islogical);
    parse(parser,varargin{:});
    params = parser.Results;
    
    valid = false;
    tries = 0;
    while ~valid
        if strcmp(distribution,'gamma')
            g_act = zeros(ge.k,1);
            for d = 1:ge.k            
                if ge.nullComponent && d == ge.k
                    g_act(d,1) = gamrnd(ge.null_shape,ge.null_scale);
                else
                    g_act(d,1) = gamrnd(ge.g_shape,ge.g_scale);
                end
            end
        elseif strcmp(distribution,'dirichlet')
             g_act = symmetricDirichlet(ge.sparsity,ge.k,1)';
        else
            error('Prior distribution %s not implemented',distribution);
        end

        valid = gestaltCheckG(g_act,ge,params.precision);
        tries = tries + 1;
        % if we cannot find a valid g, return an error code
        if tries > params.sampleRetry
            throw(MException('Gestalt:Prior:TooManyTries',sprintf('Number of tries to draw a valid g vector from the prior exceeded %d',params.sampleRetry)));            
        end
    end
end