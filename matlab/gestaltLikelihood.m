function [coefficient, exponent] = gestaltLikelihood(ge,L)
    % approximated, up to a constant
    % get L samples from a k-dimensional symmetric dirichelet prior for g
    G = symmetricDirichlet(ge.sparsity,ge.k,ge.N*L);
    coefficient = 0;
    exponent = 0;
    coeffs = zeros(L,ge.N*ge.B);
    exps = zeros(L,ge.N*ge.B);
    for s=1:L        
        for n=1:ge.N
            g = G((s-1)*ge.N+n,:)';
            Cv = componentSum(g,ge.cc);
            C = ge.obsVar * eye(ge.Dx) + ge.A * Cv * ge.A';
            for b=1:ge.B                
                x = squeeze(ge.X(n,b,:));                                
                p = mvnpdf(x,zeros(size(x)),C);
                idx = (n-1) * ge.B + b;
                [coeffs(s,idx),exps(s,idx)] = sciNot(p);                                
            end
        end 
        expsum = sum(exps(s,:),2);
        coeffprod = prod(coeffs(s,:),2);
        [co,ex] = sciNot(coeffprod);
        expsum = expsum + ex;
        if s==1
            coefficient = co;
            exponent = expsum;
        else
            [coefficient,exponent] = sumSciNot(coefficient,exponent,co,expsum);
        end
    end
    % division by L could be done
end

function [coefficient,exponent] = sciNot(a)
    [coefficient,exponent] = strread(strrep(sprintf('%E',a),'E','#'),'%f#%f');
end

function [coefficient,exponent] = sumSciNot(c1,e1,c2,e2)
    exponent = max(e1,e2);
    expdiff = abs(e1-e2);
    if(e1 > e2)
        coefficient = c1 + 10^(-expdiff) * c2;
    else
        coefficient = c2 + 10^(-expdiff) * c1;
    end
    [co,ex] = sciNot(coefficient);
    coefficient = co;
    exponent = exponent + ex;
end