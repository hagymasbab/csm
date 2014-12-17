function grad = gestaltFullLogPosteriorGrad(ge,X,V,g,z,iC)
    grad = zeros(ge.k + ge.B * ge.Dv + 1,1);
    
    % grad of g
    vv = V' * V;
    common = iC - iC * vv * iC;
    for i = 1 : ge.k
        matpart = - trace(common * ge.cc{i}) / 2;
        if ge.nullComponent && i == ge.k
            shg = ge.null_shape;
            scg = ge.null_scale;
        else
            shg = ge.g_shape;
            scg = ge.g_scale;
        end
        priorpart = (shg - 1) / g(i,1) - g(i,1) / scg;
        grad(i,1) = matpart + priorpart;
    end
    
    % grad of v
    for b = 1:ge.B
        first = ge.k + (b-1)*ge.Dv + 1;
        last = ge.k + b*ge.Dv;
        v = V(b,:)';
        x = X(b,:)';
        actgrad = - z * ( (z * ge.AA + iC) * v - ge.A' * x) / ge.obsVar;        
        grad(first:last,1) = actgrad;
    end
    
    %grad of z
    quadlin = 0;
    for b = 1:ge.B
        quadlin = quadlin + 2 * V(b,:) * ge.AA * V(b,:)' * z - X(b,:) * ge.A * V(b,:)';
    end    
    grad(end,1) = quadlin + (ge.z_shape - 1) / z - z / ge.z_scale;
end