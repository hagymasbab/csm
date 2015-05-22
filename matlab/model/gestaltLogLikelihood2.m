function ll = gestaltLogLikelihood2(ge,L,data,cholesky,varargin)
    % Mate version
    
    parser = inputParser;   
    addParameter(parser,'verbose',0,@isnumeric);  
    addParameter(parser,'loadSamples',false,@islogical);
    addParameter(parser,'randseed','leave');
    addParameter(parser,'method','intuition');
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
    siATA = ge.obsVar * inv(ATA);
    idATA = 1 / sqrt( det(ATA) );
    
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
    parfor l=1:L
        g_l = G(l,:)';
        z_l = Z(l,1);
        z2_l = Z2(l,1);
        cv = componentSum(g_l,cc);        
        if strcmp(params.method,'intuition')
            C_l = ge.obsVar * eye(ge.Dv) + z2_l * ge.A * cv * ge.A';
        elseif strcmp(params.method,'algebra')
            C_l = siATA / z2_l + cv;
            hz(l) = idATA / (z_l^ge.Dv);
        end
        % cheating
        C_l = nearestSPD(C_l);
        covariances(l,:,:) = C_l;
    end
    
    ll = 0;
    parfor n = 1:N
        ll_part = 0;
        x_n = data(n,:)';
        
        for l = 1:L                        
            if strcmp(params.method,'algebra')
                f_nl = pA * x_n / Z(l,1);
                ll_part = ll_part + hz(l) * mvnpdf(f_nl',zeros(1,ge.Dx),reshape(covariances(l,:,:),ge.Dv,ge.Dv));
            elseif strcmp(params.method,'intuition')
                ll_part = ll_part + mvnpdf(x_n',zeros(1,ge.Dx),reshape(covariances(l,:,:),ge.Dv,ge.Dv));
            end
        end
        ll = ll + log(ll_part);
    end
end