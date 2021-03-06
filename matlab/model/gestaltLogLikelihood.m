function ll = gestaltLogLikelihood(ge,L,data,varargin)

    parser = inputParser;
    addParameter(parser,'cholesky',[]);    
    addParameter(parser,'loadSamples',false,@islogical);
    addParameter(parser,'randseed','leave');      
    addParameter(parser,'verbose',0,@isnumeric);
    addParameter(parser,'scientific',true,@islogical);    
    addParameter(parser,'newway',false,@islogical);    
    parse(parser,varargin{:});        
    params = parser.Results;      

    % approximated, up to a constant
    % get L samples from the priors of g and z
    
    setrandseed(params.randseed);
    N = size(data,1);
    if ge.B ~= 1
        error('not implemented');
    end
    data = reshape(data,N,ge.Dx);
    %size(data)
    
    if params.newway
        if params.loadSamples
            % TODO check if exists
            load('prior_samples.mat'); % should contain V and Z
            if L < size(Z,1)
                Z = Z(1:L,:);
                V = V(1:L,:);
            end
        else
            V = gestaltSamplePriorV(ge,L,params.cholesky);
            Z = gamrnd(ge.z_shape,ge.z_scale,[L 1]);
        end

        sample_coeffs = zeros(L,1);
        sample_exps = zeros(L,1);
        covmat = ge.obsVar * eye(ge.Dv);
        for i=1:L
            datum_like_coeffs = zeros(N,1);
            datum_like_exps = zeros(N,1);
            act_z = Z(i,1);
            act_v = V(i,:)';
            act_mean = act_z * ge.A * act_v;
            for nn = 1:N
                [datum_like_coeffs(nn,1), datum_like_exps(nn,1)] = stableMvnpdf(data(nn,:)',act_mean,covmat,true,false);
            end
            [c,e] = prodSciNot(datum_like_coeffs',datum_like_exps');
            sample_coeffs(i,1) = c;
            sample_exps(i,1) = e;
        end

        final_coeff = 0;
        final_exp = 0;
        for i=1:L
            [final_coeff,final_exp] = sumSciNot(final_coeff,final_exp,sample_coeffs(i,1), sample_exps(i,1));
        end

        ll = scinot2log(final_coeff,final_exp);
    else
            
        if ~isempty(params.cholesky)
            for j=1:ge.k
                ge.cc{j} = params.cholesky{j}' * params.cholesky{j};                                             
            end
        end

        pA = pinv(ge.A);
        if rcond(ge.AA) < 1e-15
            iAA = pinv(ge.AA);
        else
            iAA = ge.AA \ eye(ge.Dv);
        end
        ll_coeff = 0;
        ll_exp = 0;
        ll = 0;
        if params.loadSamples
            % TODO check if exists
            load('prior_samples.mat'); % should contain G and Z
            if L < size(Z,1);
                Z = Z(1:L,:);
                G = G(1:L,:);
            end
        else
            G = gestaltSamplePriorG(ge,L);
            Z = gamrnd(ge.z_shape,ge.z_scale,[L 1]);
            save('bin/prior_samples.mat','G','Z');
        end

        log_act_l = zeros(N,1);

        parfor i=1:N
            if params.verbose > 0
                printCounter(i,'stringVal','Datapoint','maxVal',N);
            end
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
                if params.scientific
                    [pdf_coeff,pdf_exp] = stableMvnpdf(eval_site,zeros(size(eval_site)),cov_full,true,false);                
                    [z_coeff,z_exp] = sciNot(-ge.Dv * log(Z(sz,1)),true);
                    [sample_coeff,sample_exp] = prodSciNot([pdf_coeff z_coeff],[pdf_exp z_exp]);
                    [act_coeff,act_exp] = sumSciNot(act_coeff,act_exp,sample_coeff,sample_exp);
                else
                    pdfval = stableMvnpdf(eval_site,zeros(size(eval_site)),cov_full,false,false);
                    zmult = 1 / Z(sz,1)^(ge.Dv);
                    if zmult == Inf
                        sample_L = 0;
                    else
                        sample_L = pdfval * zmult;
                    end
                    act_L = act_L + sample_L;
                end
            end
            if params.scientific            
                %[ll_coeff,ll_exp] = sumSciNot(act_coeff,act_exp,ll_coeff,ll_exp);
                s = sign(act_coeff);
                act_log10 = act_exp + log10(abs(act_coeff));
                act = s * act_log10 / log10(exp(1));
                log_act_l(i,1) = act;
                %ll = ll + act;
                %ll = ll + log(act_coeff * 10^act_exp)
            else
                log_act_l(i,1) = log(act_L);
                %ll = ll + log(act_L);
            end
        end
        ll = sum(log_act_l);
        ll = ll - N * log(L);
    end
end
