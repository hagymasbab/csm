function [x,rr] = sliceSample(init,logpdf,stepsize,varargin)
    parser = inputParser;
    addParamValue(parser,'plot',false,@islogical);
    parse(parser,varargin{:});
    
    dim = size(init,1);
    pl = parser.Results.plot && dim == 1;   % only works in 1D     
    
    p_init = exp(logpdf(init));
    % sample slice height
    u = p_init * rand();
    if pl
        plot(init,u,'r-',init,u,'ro');
        pause
    end
    bounds = repmat(init,1,2);
    for d=1:dim
        % move the minimal point outside the slice
        inside = true;
        actpoint = init;
        while inside
            actpoint(d,1) = actpoint(d,1) - stepsize;
            if pl
                plot(actpoint,u,'r-',actpoint,u,'ro');
                pause
            end
            if exp(logpdf(actpoint)) < u
                bounds(d,1) = actpoint(d,1);
                inside = false;
            end
        end
        % move the maximal point outside the slice
        inside = true;
        actpoint = init;
        while inside
            actpoint(d,1) = actpoint(d,1) + stepsize;
            if pl
                plot(actpoint,u,'r-',actpoint,u,'ro');
                pause
            end
            if exp(logpdf(actpoint)) < u
                bounds(d,2) = actpoint(d,1);
                inside = false;
            end
        end
    end
    
    % sample from the approximated slice uniformly
    reject = true;
    rr = 0;
    while reject
        x = bounds(:,1) + rand(dim,1) .* (bounds(:,2)-bounds(:,1));
        if exp(logpdf(x)) >= u
            reject = false;
        else
            rr = rr + 1;
            % the approximation of the slice should be shrunk here, but I
            % don't know how to do that in more than one dimensions
            if dim == 1
                if x > init
                    bounds(1,2) = x;
                else
                    bounds(1,1) = x;
                end
            end
        end
    end