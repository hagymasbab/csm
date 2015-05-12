function grad = gestaltLogLikelihoodGradient(ge,L,data,varargin)
    parser = inputParser;   
    addParameter(parser,'verbose',0,@isnumeric);  
    addParameter(parser,'randseed','leave');      
    parse(parser,varargin{:});        
    params = parser.Results;  
    
    setrandseed(params.randseed);
    N = size(data,1);
    if ge.B ~= 1
        error('not implemented');
    end
    data = reshape(data,N,ge.Dx);
    
    G = gestaltSamplePriorG(ge,L);
    Z = gamrnd(ge.z_shape,ge.z_scale,[L 1]);
    Z_2 = Z .^ 2;
    
    pA = pinv(ge.A);
    ATA = ge.A' * ge.A;
    siATA = ge.obsVar * inv(ATA);
    idATA = 1 / sqrt( det(ATA) );
    
    M = zeros(k,ge.Dv,ge.Dv);
    for n = 1:N
        x_n = data(n,:)';
        pAx = pA * x_n;
        Li_n = 0;
        M_part = zeros(k,ge.Dv,ge.Dv);
        for l = 1:L            
            g_l = G(l,:)';
            z_l = Z(l,1);
            z_l2 = Z_2(l,1);
            
            h_l = 1 / (z_l^ge.Dv * idATA);
            f_l = pAx / z_l;
            C_v = componentSum(g_l,ge.cc);
            C_l = siATA / z_l2 + C_v;
            iC_l = inv(C_l);
            iCf = iC_l * f_l;
            
            N_f = mvnpdf(f_l',zeros(1,ge.Dv),C_l);
            scalar_term = h_f * N_f;            
            
            matrix_term = iC_l - iCf * iCf';
            
            for kk = 1:ge.k
                M_part(kk,:,:) = M_part(kk,:,:) + g_l(kk,1) * scalar_term * matrix_term;
            end
            
            Li_n = Li_n + scalar_term;            
        end
        
        M = M + M_part / Li_n;        
    end
    
    M = -M / 2;
    
    grad = cell(1,ge.k);
    for kk = 1:ge.k
        grad{kk} = zeros(ge.Dv);
        for i = 1:ge.Dv
            for j = i:ge.Dv
                grad{kk}(i,j) = 
end