function [C_m,G,Z] = msmPosteriorCovariance(x,L,randseed,loadSamples,A,B,sigma_x,sigma_v,g_shape,g_scale,z_shape,z_scale,target_acceptance)
    setrandseed(randseed);
    
    Dv = size(A,2);
    sI = eye(Dv) / sigma_v;
    sB = B / sigma_v;
    ATA = A' * A;
    sATA = ATA / sigma_x;    
    sAx = A' * x / sigma_x;
    
    if loadSamples
        load('save_msm_postcov.mat');
    else
        [G,Z] = msmHamiltonianGZ(x,L,50,'leave',A,B,sigma_x,sigma_v,g_shape,g_scale,z_shape,z_scale,target_acceptance);
%         G = real(G);
%         Z = real(Z);
        save('bin/save_msm_postcov.mat','G','Z');
    end
    
    Cml = zeros(L,Dv,Dv);        
    mml = zeros(L,Dv);

    for i = 1:L      
        z2ATA = Z(i)^2 * sATA;        
        nCm = stableInverse(z2ATA + sI);
        %[~,e] = cholcov(nCm)
        %Cm = nearestSPD(nCm);        
        Cm = nCm;
        mm = Cm * (Z(i)*sAx + sB*G(i,:)');        
        Cml(i,:,:) = Cm;        
        mml(i,:) = mm';
    end
    
    C_m = squeeze(mean(Cml,1)) + cov(mml,1);    

end