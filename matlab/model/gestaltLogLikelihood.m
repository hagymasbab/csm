function ll = gestaltLogLikelihood(ge,L,data,cholesky)
    % approximated, up to a constant
    % get L samples from the priors of g and z
    
    if ~isempty(cholesky)
        for j=1:ge.k
            ge.cc{j} = cholesky{j}' * cholesky{j};                                             
        end
    end
    
    pA = pinv(ge.A);
    siAA = ge.obsVar * inv(ge.A*ge.A');
    ll = 0;

    G = gestaltSamplePriorG(ge,L);
    Z = gamrnd(ge.z_shape,ge.z_scale,[L 1]);
    for i=1:size(data,1)
        % TODO B > 1
        x = data(i,:)';
        Ax = pA * x;
        act_L = 0;
        for sz=1:L
            cov_left = siAA / Z(sz,1)^2;
            eval_site = Ax / Z(sz,1);
            act_LG = 0;
            for sg=1:L
                g = G(sg,:)';
                cov_full = cov_left + componentSum(g,ge.cc);
                act_LG = act_LG + normpdf(eval_site',zeros(size(eval_site)),cov_full);
            end
            act_L = act_L + act_LG / Z(sz,1)^(ge.Dv);
        end
        ll = ll + log(act_L);
    end
%             batch_coeffs = zeros(1,ge.B);
%             batch_exps = zeros(1,ge.B);
%             g = G((i-1)*L+s,:)';
%             Cv = componentSum(g,ge.cc);
%             %C = ge.obsVar * eye(ge.Dx) + ge.A * Cv * ge.A';
%             C = ge.obsVar * iAA + Cv;
%             [~,err] = chol(C);
%             if err == 0 && isequal(C,C')                                          
%                 for b=1:ge.B                
%                     x = squeeze(ge.X(n,b,:)); 
%                     x = pA * x;
%                     p = mvnpdf(x,zeros(size(x)),C);                               
%                     [batch_coeffs(1,b),batch_exps(1,b)] = sciNot(p);                                
%                 end            
%             end
%             [samp_coeffs(1,s),samp_exps(1,s)] = sciProd(batch_coeffs,batch_exps);
%         end 
%         [datum_coeff,datum_exp] = sciSum(samp_coeffs,samp_exps);
%         ll = ll + log10(datum_coeff) + datum_exp;       
%     end
end
