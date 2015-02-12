function V = gestaltSamplePriorV(ge,L,cholesky)
    if ge.B ~= 1
        error('not implemented');
    end
    if ~isempty(cholesky)
        for j=1:ge.k
            ge.cc{j} = cholesky{j}' * cholesky{j};                                             
        end
    end
    G = gestaltSamplePriorG(ge,L);    
    V = zeros(L,ge.Dv);
    for i = 1:L
        Cv = componentSum(G(i,:)',ge.cc);
        V(i,:) = mvnrnd(zeros(1,ge.Dv),Cv);
    end
end