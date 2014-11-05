function LM = signlog(M)
    LM = M;
    LM(M>0) = log(M(M>0));
    LM(M<0) = -log(-M(M<0));
end