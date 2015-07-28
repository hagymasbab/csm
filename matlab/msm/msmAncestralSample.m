function [X,V] = msmAncestralSample(g,z,randseed,A,B,sigma_x,sigma_v)
    setrandseed(randseed);
    Dv = size(A,2);
    Dx = size(A,1);
    V = mvnrnd((B * g)',sigma_v*eye(Dv))';        
    X = mvnrnd( (z*A * V)',sigma_x*eye(Dx));
end