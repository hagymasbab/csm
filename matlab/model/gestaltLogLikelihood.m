function ll = gestaltLogLikelihood(ge,L,data,cholesky)
    % approximated, up to a constant
    % get L samples from the priors of g and z
    
    if ~isempty(cholesky)
        for j=1:ge.k
            ge.cc{j} = cholesky{j}' * cholesky{j};                                             
        end
    end
    
    pA = pinv(ge.A);
    iAA = ge.AA \ eye(ge.Dv);
    ll = 0;

    G = gestaltSamplePriorG(ge,L);
    Z = gamrnd(ge.z_shape,ge.z_scale,[L 1]);
    for i=1:size(data,1)
        % TODO B > 1
        x = data(i,:)';
        Ax = pA * x;
        act_L = 0;
        for sz=1:L
            %ge.obsVar/Z(sz,1)^2
            cov_left = (ge.obsVar/Z(sz,1)^2) * iAA;
            %[~,e] = chol(cov_left);
            %fprintf('siAA/z^2 chol: %d det %E det_stab %E\n',e,det(cov_left),stableLogdet(cov_left));                        
            eval_site = Ax / Z(sz,1);
            act_LG = 0;
            %for sg=1:L
                sg = sz;
                g = G(sg,:)';
                cv = componentSum(g,ge.cc);
                %[~,e] = chol(cv);
                %fprintf('Cv chol: %d det %E det_stab %E\n',e,det(cv),stableLogdet(cv));    
                cov_full = cov_left + cv;
                %cov_full = cv;
                %[~,e] = chol(cov_full);
                %fprintf('Full chol: %d det %E det_stab %E\n',e,det(cov_full),stableLogdet(cov_full));    
                %pdfval = mvnpdf(eval_site',zeros(size(eval_site))',cov_full);
                pdfval = stableMvnpdf(eval_site,zeros(size(eval_site)),cov_full);
                %pause
                act_LG = act_LG + pdfval;
            %end
            act_L = act_L + act_LG / Z(sz,1)^(ge.Dv);
        end
        if act_L ~= 0
            ll = ll + log(act_L);
        end
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
