function ll = gestaltLogLikelihood(ge,L,data,cholesky,loadSamples,randseed,scientific)
    % approximated, up to a constant
    % get L samples from the priors of g and z
    
    setrandseed(randseed);
    N = size(data,1);
    % TODO nem megy B~=1-re
    data = reshape(data,N,ge.Dx);
        
    if ~isempty(cholesky)
        for j=1:ge.k
            ge.cc{j} = cholesky{j}' * cholesky{j};                                             
        end
    end
    
    pA = pinv(ge.A);
    iAA = ge.AA \ eye(ge.Dv);
    %iAA = inv(ge.AA);
    ll_coeff = 0;
    ll_exp = 0;
    ll = 0;
    if loadSamples
        % TODO check if exists
        load('prior_samples.mat'); % should contain G and Z
        if L < size(Z,1);
            Z = Z(1:L,:);
            G = G(1:L,:);
        end
    else
        G = gestaltSamplePriorG(ge,L);
        Z = gamrnd(ge.z_shape,ge.z_scale,[L 1]);
    end
    
    for i=1:N
        % TODO B > 1
        x = data(i,:)';
        Ax = pA * x;
        act_L = 0;
        act_coeff = 0;
        act_exp = 0;
        for sz=1:L
            cov_left = (ge.obsVar/Z(sz,1)^2) * iAA;                                
            eval_site = Ax / Z(sz,1);
            g = G(sz,:)';
            cv = componentSum(g,ge.cc);           
            cov_full = cov_left + cv;            
            if scientific
                [pdf_coeff,pdf_exp] = stableMvnpdf(eval_site,zeros(size(eval_site)),cov_full,true);
                [z_coeff,z_exp] = sciNot(-ge.Dv * log(Z(sz,1)),true);
                [sample_coeff,sample_exp] = prodSciNot([pdf_coeff z_coeff],[pdf_exp z_exp]);
                [act_coeff,act_exp] = sumSciNot(act_coeff,act_exp,sample_coeff,sample_exp);
            else
                pdfval = stableMvnpdf(eval_site,zeros(size(eval_site)),cov_full,false);
                zmult = 1 / Z(sz,1)^(ge.Dv);
                sample_L = pdfval * zmult;
                act_L = act_L + sample_L;
            end
        end
        if scientific            
            %[ll_coeff,ll_exp] = sumSciNot(act_coeff,act_exp,ll_coeff,ll_exp);
            s = sign(act_coeff);
            act_log10 = act_exp + log10(abs(act_coeff));
            act = s * act_log10 / log10(exp(1));
            ll = ll + act;
            %ll = ll + log(act_coeff * 10^act_exp)
        else
            ll = ll + log(act_L);
        end
    end
    ll = - ll * N * log(L);
end
