function grad = gestaltPostGGrad(g,v,ge)
    dimred = false;
    if size(g,1) < ge.k
        g = [g;1-sum(g)];
        dimred = true;
    end
    vv = v*v';
    Cv = componentSum(g,ge.cc);
    leftmat = Cv \ (eye(ge.Dv) - vv / Cv);
    grad = zeros(ge.k,1);
    for i=1:ge.k             
        rightmat = ge.cc{i};
        grad(i,1) = (-1/2) * reshape(leftmat,1,ge.Dv*ge.Dv) * reshape(rightmat,ge.Dv*ge.Dv,1) + (ge.sparsity - 1)/g(i,1);
    end
    
    if dimred
        grad = grad(1:size(grad,1)-1,:);
    end
end
        