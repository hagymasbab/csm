function grad = gestaltLogPostGradient(g,v,x,ge,gOff)
    grad = zeros(ge.k+ge.Dv,1);
    % check whether we already have the inverse covariance for the given g,
    % and store it if not
%     if ge.lastG == g
%         Cv = ge.Cv;
%         fprintf('e');
%     else
%         ge.lastG = g;
%         Cv = componentSum(g,ge.cc);
%         ge.Cv = Cv;
%     end
    Cv = componentSum(g,ge.cc);
    % v
    %grad(ge.k+1:size(grad,1),1) = -ge.gradMat * (ge.A*v - x ) - iCv * v;
    grad(ge.k+1:size(grad,1),1) = -ge.gradMat * (ge.A*v - x ) - Cv \ v;
        
    % g
    if ~gOff
        vv = v'*v;
        leftmat = iCv * (eye(ge.Dx) - vv*iCv);
        for i=1:ge.k
            if g(i,1) > 1
                grad(i,1) = -Inf;
            elseif g(i,1) < 0
                grad(i,1) = -Inf;
            else         
                rightmat = ge.cc{i};
                grad(i,1) = (-1/2) * reshape(leftmat,1,ge.Dv*ge.Dv) * reshape(rightmat,ge.Dv*ge.Dv,1) + (ge.sparsity - 1)/g(i,1);
            end
        end        
    end
end