function [x,rr] = sliceSample(init,logpdf,stepsize,varargin)
    parser = inputParser;
    addParamValue(parser,'plot',false,@islogical);
    addParamValue(parser,'verbose',false,@islogical);
    addParamValue(parser,'sampleRetry',100,@isnumeric);
    parse(parser,varargin{:});
    verb = parser.Results.verbose;
    retry = parser.Results.sampleRetry;
    
    dim = size(init,1);
    pl = parser.Results.plot && dim == 1;   % only works in 1D     
    
    logp_init = logpdf(init);
    %p_init = exp(logpdf(init));
    % sample slice height
    %u = p_init * rand();
    log_u = logp_init + log(rand());
    if pl
        %plot(init,u,'r-',init,u,'ro');
        plot(init,log_u,'r-',init,log_u,'ro');
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
                %plot(actpoint,u,'r-',actpoint,u,'ro');
                plot(actpoint,log_u,'r-',actpoint,log_u,'ro');
                pause
            end
            %if exp(logpdf(actpoint)) <= u
            if logpdf(actpoint) <= log_u
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
                %plot(actpoint,u,'r-',actpoint,u,'ro');
                plot(actpoint,log_u,'r-',actpoint,log_u,'ro');
                pause
            end
            %if exp(logpdf(actpoint)) <= u
            if logpdf(actpoint) <= log_u
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
        %if exp(logpdf(x)) >= u
        if logpdf(x) >= log_u
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
        if rr > retry
           rr = -1;
           return;
        end
    end
end