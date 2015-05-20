function [C_c,C_m] = posteriorCovariances(x,A,sigma_x,sigma_v,shape_g,scale_g,shape_z,scale_z,cc,B,L,randseed)
    setrandseed(randseed);
        
    k = size(B,2);
    Dx = size(A,1);
    Dv = size(A,2);
    ATA = A' * A;
    sATA = ATA / sigma_x;    
    sAx = A' * x / sigma_x;
    sI = eye(Dv) / sigma_v;
    sB = B / sigma_v;
    
    G = gamrnd(shape_g,scale_g,[L k]);
    Z = gamrnd(shape_z,scale_z,[L 1]);
    
    Ccl = zeros(L,Dv,Dv);
    Cml = zeros(L,Dv,Dv);    
    mcl = zeros(L,Dv);
    mml = zeros(L,Dv);
    
    for i = 1:L
        z2ATA = Z(i)^2 * sATA;
        Cvl = componentSum(G(i,:)',cc);
        iCvl = inv(Cvl);
        Cc = inv(z2ATA + iCvl);
        mc = Z(i) * Cc * sAx;
        Cm = inv(z2ATA + sI);
        mm = Cm * (Z(i)*sAx + sB*G(i,:)');
        Ccl(i,:,:) = Cc;
        Cml(i,:,:) = Cm;
        mcl(i,:) = mc';
        mml(i,:) = mm';
    end
    
    C_c = squeeze(mean(Ccl,1)) + cov(mcl,1);
    C_m = squeeze(mean(Cml,1)) + cov(mml,1);
end