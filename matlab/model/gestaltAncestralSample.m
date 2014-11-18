function [X,V] = gestaltAncestralSample(ge,g,z,precision,positive)
    if ~precision
        Cv = componentSum(g,ge.cc);
    else
        Cv = inv(componentSum(g,ge.pc));
    end
    
    V = mvnrnd(zeros(ge.B,ge.Dv),Cv);
    if positive
        V = abs(V);        
    end

    means = reshape(V,ge.B,ge.Dv);
    means = z * means * ge.A';
    X = mvnrnd(means,ge.obsVar*eye(ge.Dx));
end