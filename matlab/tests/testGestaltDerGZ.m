function testGestaltDerGZ(formula,randseed)
        
    function zg = zgrad(g,z,x,ge)
        gr = gestaltLogPostGZGrad(g,z,x,ge);
        zg = gr(end);
    end

    function gg = ggrad(g,z,x,ge)
        gr = gestaltLogPostGZGrad(g,z,x,ge);
        gg = gr(1:ge.k);
    end

    function zg = msm_zgrad(g,z,x,A,B,sigma_x,sigma_v,g_shape,g_scale,z_shape,z_scale)
        gr = msmLogPostGZGrad(g,z,x,A,B,sigma_x,sigma_v,g_shape,g_scale,z_shape,z_scale);
        zg = gr(end);
    end

    function gg = msm_ggrad(g,z,x,A,B,sigma_x,sigma_v,g_shape,g_scale,z_shape,z_scale)
        gr = msmLogPostGZGrad(g,z,x,A,B,sigma_x,sigma_v,g_shape,g_scale,z_shape,z_scale);
        gg = gr(1:end-1);
    end

    function Cx = postcov(g1,g,z,ge,kk)
        g(kk) = g1;
        Cv = componentSum(g,ge.cc);
        Cx = ge.obsVar * eye(ge.Dx) + z^2 * ge.A * Cv * ge.A';
        Cx = Cx(1,1);
    end

    function gr = postcov_derg(g1,g,z,ge,kk)
        g(kk) = g1;
        gr = z^2 * ge.A * ge.cc{kk} * ge.A';
    end

    function gr = postcov_derz(g,z,ge)
        Cv = componentSum(g,ge.cc);
        gr = 2 * z * ge.A * Cv * ge.A';
        gr = gr(1,1);
    end

    setrandseed(randseed);
    ge = gestaltCreate('temp','Dx',64,'k',2,'N',1,'filters','OF','obsVar',0.5,'g_shape',1,'g_scale',1,'z_shape',2,'z_scale',2, ...
        'generateComponents',true,'generateData',true);    
    x = reshape(ge.X(1,1,:),ge.Dx,1);
    B = randn(ge.Dv,ge.k);
    sigma_v = 0.5;
    
    if strcmp(formula,'post')
        a = @(h) gestaltLogPostGZ(h(1:ge.k),h(end),x,ge);
        b = @(h) gestaltLogPostGZGrad(h(1:ge.k),h(end),x,ge);
        init = ones(ge.k+1,1);
    elseif strcmp(formula,'zprior')
        a = @(z) (ge.z_shape - 1) * log(z) - z / ge.z_scale;
        b = @(z) (ge.z_shape-1)/z - 1/ge.z_scale;
        init = 1;
    elseif strcmp(formula,'zpost')
        a = @(z) gestaltLogPostGZ(ones(ge.k,1),z,x,ge);
        b = @(z) zgrad(ones(ge.k,1),z,x,ge);
        init = 1;
    elseif strcmp(formula,'gpost')
        a = @(g) gestaltLogPostGZ(g,1,x,ge);
        b = @(g) ggrad(g,1,x,ge);
        init = ones(ge.k,1);
    elseif strcmp(formula,'cxg')
        a = @(g) postcov(g,ones(ge.k,1),1,ge,1);
        b = @(g) postcov_derg(g,ones(ge.k,1),1,ge,1);
        init = 1;
    elseif strcmp(formula,'cxz')
        a = @(z) postcov(1,ones(ge.k,1),z,ge,1);
        b = @(z) postcov_derz(ones(ge.k,1),z,ge);
        init = 1;
    elseif strcmp(formula,'msm')        
        a = @(h) msmLogPostGZ(h(1:ge.k),h(end),x,ge.A,B,ge.obsVar,sigma_v,ge.g_shape,ge.g_scale,ge.z_shape,ge.z_scale);
        b = @(h) msmLogPostGZGrad(h(1:ge.k),h(end),x,ge.A,B,ge.obsVar,sigma_v,ge.g_shape,ge.g_scale,ge.z_shape,ge.z_scale);
        init = ones(ge.k+1,1);
    elseif strcmp(formula,'msmzpost')
        a = @(z) msmLogPostGZ(ones(ge.k,1),z,x,ge.A,B,ge.obsVar,sigma_v,ge.g_shape,ge.g_scale,ge.z_shape,ge.z_scale);  
        b = @(z) msm_zgrad(ones(ge.k,1),z,x,ge.A,B,ge.obsVar,sigma_v,ge.g_shape,ge.g_scale,ge.z_shape,ge.z_scale);
        init = 1;
    elseif strcmp(formula,'msmgpost')
        a = @(g) msmLogPostGZ(g,1,x,ge.A,B,ge.obsVar,sigma_v,ge.g_shape,ge.g_scale,ge.z_shape,ge.z_scale);  
        b = @(g) msm_ggrad(g,1,x,ge.A,B,ge.obsVar,sigma_v,ge.g_shape,ge.g_scale,ge.z_shape,ge.z_scale);
        init = ones(ge.k,1);
    end
    
    disc = checkDerivative(a,b,init,false)
end