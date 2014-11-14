function [cc,coeffs] = templates2covariances(templates,A)
    %iA = pinv(A);
    iA = A';
    k = size(templates,2);
    Dx = size(templates{1},1)^2;
    cc = cell(1,k);
    coeffs = zeros(k,Dx);
    for t = 1:k
        act_coeff = iA * templates{t}(:);
        mean_coeff = mean(act_coeff);
        %threshold = abs(mean_coeff);
        threshold = (max(abs(act_coeff)) + abs(mean_coeff)) / 2;
        c_act = zeros(Dx);
        for i = 1:Dx            
            for j = i+1:Dx
%                 c_act(i,j) = abs(abs(act_coeff(i,1)) - abs(act_coeff(j,1)));
                
                if (abs(act_coeff(i,1)) > threshold) && (abs(act_coeff(j,1)) > threshold)
                    c_act(i,j) = 1;
                end
                
                if act_coeff(i,1)*act_coeff(j,1) < 0
                    c_act(i,j) = -c_act(i,j);
                end
                
                c_act(j,i) = c_act(i,j);
            end
            c_act(i,i) = sum(abs(c_act(i,:))) + 0.1;
        end
        cc{t} = c_act;
        coeffs(t,:) = act_coeff';
    end
end