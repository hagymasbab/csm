function ll = gestaltLogLikelihood2(ge,L,data,cholesky,varargin)    
    parser = inputParser;   
    addParameter(parser,'verbose',0,@isnumeric);  
    addParameter(parser,'loadSamples',false,@islogical);
    addParameter(parser,'randseed','leave');
    addParameter(parser,'method','intuition');
    addParameter(parser,'scientific','false');
    parse(parser,varargin{:});        
    params = parser.Results;  
    
    setrandseed(params.randseed);
    N = size(data,1);
    if ge.B ~= 1 || ge.nullComponent
        error('not implemented');
    end
    data = reshape(data,N,ge.Dx);
    
    cc = cell(1,ge.k);    
    for i=1:ge.k
        cc{i} = cholesky{i}' * cholesky{i};
    end  
    pA = pinv(ge.A);
    ATA = ge.A' * ge.A;
    siATA = ge.obsVar * stableInverse(ATA);
    
    if strcmp(params.method,'algebra')
        logConstant = -N * (log(L) + stableLogdet(ATA) / 2);
    elseif strcmp(params.method,'intuition')
        logConstant = -N * log(L);
    end
    
    if params.loadSamples
        load('prior_samples.mat');
    else
        G = gestaltSamplePriorG(ge,L);
        Z = gamrnd(ge.z_shape,ge.z_scale,[L 1]);
        save('bin/prior_samples.mat','G','Z');
    end
    Z2 = Z .^ 2;
    
    covariances = zeros(L,ge.Dv,ge.Dv);
    hz = zeros(L,1);
    loghz = zeros(L,1);
    for l=1:L
        g_l = G(l,:)';
        z_l = Z(l,1);
        z2_l = Z2(l,1);
        cv = componentSum(g_l,cc);    
        %issymmetric(cv)
        if strcmp(params.method,'intuition')
            c_right = z2_l * (ge.A * cv * ge.A');
            c_right = (c_right + c_right') / 2;
            %issymmetric(c_right)
            Cl = ge.obsVar * eye(ge.Dv) + c_right;
        elseif strcmp(params.method,'algebra')
            error('algebraic approach not implemented');
            Cl = siATA / z2_l + cv;
            hz(l) = 1 / (z_l^ge.Dv);
        end
        % cheating
%         [~,e] = cholcov(c_right)    
%         issymmetric(Cl)
        nCl = nearestSPD(Cl);
%         [~,e] = cholcov(c_right)            
%         issymmetric(nCl)
        covariances(l,:,:) = nCl;
    end
    
    ll = 0;
    for n = 1:N
        ll_part = 0;
        ll_part_coeff = 0;
        ll_part_expo = 0;
        x_n = data(n,:)';
        
        for l = 1:L                        
            if strcmp(params.method,'algebra')
                f_nl = pA * x_n / Z(l,1);
                ll_part = ll_part + hz(l) * mvnpdf(f_nl',zeros(1,ge.Dx),reshape(covariances(l,:,:),ge.Dv,ge.Dv));
            elseif strcmp(params.method,'intuition')
                [ll_act_coeff,ll_act_expo] = stableMvnpdf(x_n,zeros(ge.Dx,1),reshape(covariances(l,:,:),ge.Dv,ge.Dv),true,false);
                [ll_part_coeff,ll_part_expo] = sumSciNot(ll_part_coeff,ll_part_expo,ll_act_coeff,ll_act_expo);
            end
        end
        %ll = ll + log(ll_part);
        ll = ll + scinot2log(ll_part_coeff,ll_part_expo);
    end
    
    ll = ll + logConstant;
end