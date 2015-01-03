function grad = gestaltFullLogPosteriorGrad(ge,X,V,g,z,iC)
    if ge.B ~= 1
        error('Gradient is only valid for B=1');
    end
    
    if isempty(iC)
        Cv = sparse(componentSum(g,ge.cc));
        iC = inv(Cv);
    end
    
    grad = zeros(ge.k + ge.B * ge.Dv + 1,1);
    
    % grad of g
    vv = V' * V;
    common = ge.B * iC - iC * vv * iC;
    for i = 1 : ge.k
        matpart = - trace(common * ge.cc{i}) / 2;
        if ge.nullComponent && i == ge.k
            shg = ge.null_shape;
            scg = ge.null_scale;
        else
            shg = ge.g_shape;
            scg = ge.g_scale;
        end
        priorpart = (shg - 1) / g(i,1) - 1 / scg;
        grad(i,1) = matpart + priorpart;
    end
    
    % grad of v
    for b = 1:ge.B
        first = ge.k + (b-1)*ge.Dv + 1;
        last = ge.k + b*ge.Dv;
        v = V(b,:)';
        x = X(b,:)';
        constant_term = (z / ge.obsVar) * ge.A' * x;
        linear_term = ((z^2 / ge.obsVar) * ge.AA + iC) * v;            
        grad(first:last,1) = constant_term - linear_term;
    end
    
    %grad of z
    quadlin = 0;
    for b = 1:ge.B
        quadlin = quadlin + X(b,:) * ge.A * V(b,:)' - V(b,:) * ge.AA * V(b,:)' * z;
    end    
    grad(end,1) = (1/ge.obsVar) * quadlin + (ge.z_shape - 1) / z - 1 / ge.z_scale;
end