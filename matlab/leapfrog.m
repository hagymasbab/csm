function [p,q,ptraj,qtraj,refl] = leapfrog(p,q,grad,stepSize,lfSteps,bounds)
    dim = size(q,1);
    refl = 0;
    boundnum = size(bounds,1);
    qtraj = zeros(dim,lfSteps+1);
    ptraj = zeros(dim,lfSteps+1);
    qtraj(:,1) = q;
    ptraj(:,1) = p;
    p = p - stepSize * grad(q)/2;
    for l=1:lfSteps
        q = q + stepSize * p;           
        if l<lfSteps
            p = p - stepSize * grad(q);
        else
            p = p - stepSize * grad(q)/2;
        end
        
        for b=1:boundnum
            bdim = bounds(b,1);
            within = false;
            while(~within)
                if q(bdim) < bounds(b,2)
                    q(bdim) = bounds(b,2) + (bounds(b,2) - q(bdim));
                    p(bdim) = -p(bdim);
                    refl = refl + 1;
                elseif q(bdim) > bounds(b,3)
                    q(bdim) = bounds(b,3) - (q(bdim) - bounds(b,3));
                    p(bdim) = -p(bdim);
                    refl = refl + 1;
                else
                    within = true;
                end
            end                
        end
        
        qtraj(:,l+1) = q;
        ptraj(:,l+1) = p;
    end   
end