function [cc,coeffs] = templates2covariances(templates,A)
    iA = pinv(A);
    k = size(templates,2);
    Dx = size(templates{1},1)^2;
    cc = cell(1,k);
    coeffs = zeros(k,Dx);
    for t = 1:k
        act_coeff = iA * templates{t}(:);
        c_act = zeros(Dx);
        for i = 1:Dx            
            for j = i+1:Dx
                c_act(i,j) = min(act_coeff(i,1),act_coeff(j,1));
                c_act(j,i) = c_act(i,j);
            end
            %c_act(i,i) = sum(c_act(i,:)) + 0.1;
        end
        cc{t} = c_act;
        coeffs(t,:) = act_coeff';
    end
end