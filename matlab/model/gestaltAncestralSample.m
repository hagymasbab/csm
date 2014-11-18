function [X,V] = gestaltAncestralSample(ge,g,z,precision)
    if ~precision
        Cv = componentSum(g,ge.cc);
    else
        Cv = inv(componentSum(g,ge.pc));
    end
    %V = mvnrnd(zeros(ge.B,ge.Dv),Cv);                                
    V = abs(mvnrnd(zeros(ge.B,ge.Dv),Cv));

    means = reshape(V,ge.B,ge.Dv);
    means = z * means * ge.A';
    X = mvnrnd(means,ge.obsVar*eye(ge.Dx));
end