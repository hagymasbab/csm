function V = gestaltSamplePriorV(ge,L)
    if ge.B ~= 1
        error('not implemented');
    end
    G = gestaltSamplePriorG(ge,L);    
    V = zeros(L,ge.Dv);
    for i = 1:L
        Cv = componentSum(G(i,:)',ge.cc);
        V(i,:) = mvnrnd(zeros(1,ge.Dv),Cv);
    end
end