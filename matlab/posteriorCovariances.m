function C_c = posteriorCovariances(x,ge,L,randseed,loadSamples)
    %function [C_c,C_m] = posteriorCovariances(x,A,sigma_x,sigma_v,shape_g,scale_g,shape_z,scale_z,cc,B,L,randseed)
    
    setrandseed(randseed);
    
    A = ge.A;
    Dv = size(A,2);
    ATA = A' * A;
    sATA = ATA / ge.obsVar;    
    sAx = A' * x / ge.obsVar;
%     sI = eye(Dv) / sigma_v;
%     sB = B / sigma_v;
    
    % take posterior samples
    % TODO load from multiple files
    if loadSamples
        load('save_postcov.mat');
    else
        [G,Z] = gestaltHamiltonianGZ(x,ge,L,50,'leave');
        save('bin/save_postcov.mat','G','Z');
    end
    
    Ccl = zeros(L,Dv,Dv);
    mcl = zeros(L,Dv);
    
%     Cml = zeros(L,Dv,Dv);        
%     mml = zeros(L,Dv);
    
    for i = 1:L
        z2ATA = Z(i)^2 * sATA;
        Cvl = componentSum(G(i,:)',ge.cc);
        iCvl = stableInverse(Cvl);
        Cc = stableInverse(z2ATA + iCvl);
        mc = Z(i) * Cc * sAx;
        Ccl(i,:,:) = Cc;
        mcl(i,:) = mc';
        
%         Cm = stableInverse(z2ATA + sI);
%         mm = Cm * (Z(i)*sAx + sB*G(i,:)');        
%         Cml(i,:,:) = Cm;        
%         mml(i,:) = mm';
    end
    
    C_c = squeeze(mean(Ccl,1)) + cov(mcl,1);
%     C_m = squeeze(mean(Cml,1)) + cov(mml,1);    
end