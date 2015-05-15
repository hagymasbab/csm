function ll = gestaltLogLikelihood2(ge,L,data,cholesky,varargin)
    % Mate version
    
    parser = inputParser;   
    addParameter(parser,'verbose',0,@isnumeric);  
    addParameter(parser,'loadSamples',false,@islogical);
    addParameter(parser,'randseed','leave');      
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
    
    if params.loadSamples
        load('prior_samples.mat');
    else
        G = gestaltSamplePriorG(ge,L);
        Z = gamrnd(ge.z_shape,ge.z_scale,[L 1]);
        save('bin/prior_samples.mat','G','Z');
    end
    Z2 = Z .^ 2;
    
    covariances = zeros(L,ge.Dv,ge.Dv);
    for l=1:L
        g_l = G(l,:)';
        z2_l = Z2(l,1);
        cv = componentSum(g_l,cc);
        covariances(l,:,:) = ge.obsVar * eye(ge.Dv) + z2_l * ge.A * cv * ge.A';
    end
    
    ll = 0;
    for n = 1:N
        ll_part = 0;
        x_n = data(n,:)';
        for l = 1:L
            ll_part = ll_part + mvnpdf(x_n',zeros(1,ge.Dx),reshape(covariances(l,:,:),ge.Dv,ge.Dv));
        end
        ll = ll + log(ll_part);
    end
end