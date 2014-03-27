classdef GestaltLogPostGrad
    properties (SetAccess = private, GetAccess = private)
        Cv
        lastG
        gradMat
    end
    
    methods
        function obj = GestaltLogPostGrad(ge)
            obj.Cv = NaN * ones(ge.Dv,ge.Dv);
            obj.gradMat = -ge.A' * (1/ge.obsVar) * eye(ge.Dx);
            obj.lastG = NaN * ones(ge.k,1);
        end
        
        function grad = statefulGrad(obj,g,v,x,ge)
            if obj.lastG ~= g
                obj.lastG = g;
                obj.Cv = componentSum(g,ge.cc);
            end
            grad = zeros(ge.k+ge.Dv,1);
            grad(ge.k+1:ge.Dv+ge.k,1) = obj.gradMat * (ge.A*v - x ) - obj.Cv \ v;
            
            % WARNING now it always gives zeros for the g dimensions for
            % the sake of speed, if you need thos too, this part should be
            % included
%             if ~gOff
%                 vv = v'*v;
%                 leftmat = iCv * (eye(ge.Dx) - vv*iCv);
%                 for i=1:ge.k
%                     if g(i,1) > 1
%                         grad(i,1) = -Inf;
%                     elseif g(i,1) < 0
%                         grad(i,1) = -Inf;
%                     else         
%                         rightmat = ge.cc{i};
%                         grad(i,1) = (-1/2) * reshape(leftmat,1,ge.Dv*ge.Dv) * reshape(rightmat,ge.Dv*ge.Dv,1) + (ge.sparsity - 1)/g(i,1);
%                     end
%                 end
%             end
        end
    end
end