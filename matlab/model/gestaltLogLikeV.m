function loglike = gestaltLogLikeV(V,g,ge,precision,iC)
    if ndims(V) == 3
        V = reshape(V,ge.B,ge.Dv);
    end
    B = size(V,1);        
    
    if ~precision
        Cv = componentSum(g,ge.cc);
        Cv = sparse(Cv);
        [~,err] = chol(Cv);
        if err > 0
            loglike = -Inf;
            return;
        end
    else
        P = componentSum(g,ge.pc);
    end        
    
    quad = 0;
    for b=1:B
        vb = V(b,:)';
        
        if ~precision
            if isempty(iC)
                quad = quad + vb' * (Cv \ vb);
            else
                quad = quad + vb' * iC * vb;
            end
            %this should be faster
%             opts.LT = true;
%             opts.UT = false;
%             temp = linsolve(U',vb,opts);
%             opts.LT = false;
%             opts.UT = true;
%             rightvec = linsolve(U,temp,opts);
%             quad = quad + vb' * rightvec;
        else
            quad = quad + vb' * P * vb;
        end
    end
    if ~precision        
        %logdet = log(det(Cv)); % this can be numerically very unstable
        logdet = sum(log(eig(Cv)));
        loglike = (-1/2) * ( B* logdet + quad );
    else
        loglike = (-1/2) * ( B* log(1/det(P)) + quad );
    end
end