function [x,rr] = sliceSample(init,logpdf,stepsize)
    dim = size(init,1);
    p_init = exp(logpdf(init));
    % sample slice height
    u = p_init * rand();
    bounds = repmat(init,1,2);
    for d=1:dim
        % move the minimal point outside the slice
        inside = true;
        actpoint = init;
        while inside
            actpoint(d,1) = actpoint(d,1) - stepsize;
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
        end
    end