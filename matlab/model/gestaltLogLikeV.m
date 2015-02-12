function loglike = gestaltLogLikeV(V,g,ge,precision,iC,Cv)
    if ndims(V) == 3
        V = reshape(V,ge.B,ge.Dv);
    end
    B = size(V,1);        
    
    illConditioned = false;
    if ~precision
        if isempty(Cv)
            Cv = componentSum(g,ge.cc);            
        end
        if rcond(Cv) < 1e-15
            illConditioned = true;
        end
        if ge.sparseComponents
            Cv = sparse(Cv);
        end
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
                if ~illConditioned
                    quad = quad + vb' * (Cv \ vb);
                else
                    quad = quad + vb' * pinv(Cv) * vb;
                end
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
        logDet = stableLogdet(Cv);
        loglike = (-1/2) * ( B* logDet + quad );
    else
        loglike = (-1/2) * ( B* log(1/det(P)) + quad );
    end
end