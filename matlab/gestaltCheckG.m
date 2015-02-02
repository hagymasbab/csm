function good = gestaltCheckG(g,ge,precision)
    good = true;
    if ~precision
        CvP = componentSum(g,ge.cc);
    else
        CvP = componentSum(g,ge.pc);
    end
    if rcond(CvP) < 1e-15
        good = false;
        return;
    end
    if ~precision
        postP = (1/ge.obsVar) * ge.AA + inv(CvP);                
    else
        postP = (1/ge.obsVar) * ge.AA + CvP;                
    end
    if rcond(postP) < 1e-15
        good = false;
        return;
    end
    
    %postC = inv(postP);
    [~,err] = cholcov(postP);
    %if det(CvP) > 0 && det(postC) > 0 && err == 0                
    if err ~= 0                
        good = false;                
    end
end