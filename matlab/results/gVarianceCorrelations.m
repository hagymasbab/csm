function gVarianceCorrelations(g_scales,move_mean)
    %cc = eyes(2,64,1);
    cc = randomCovariances(2,64,'noiseLevel',10);
    ge = gestaltCreate('temp','Dx',64,'k',2,'filters','gabor_4or','obsVar',0.7,'g_shape',1,'g_scale',0.1,'cc',cc, ...
        'z_shape',2,'z_scale',2,'N',1,'generateComponents',false,'generateData',true);

%     ge = gestaltCreate('temp','Dx',64,'k',2,'filters','gabor_4or','obsVar',0.7,'g_shape',1,'g_scale',0.1, ...
%         'z_shape',2,'z_scale',2,'N',1000,'generateComponents',true,'generateData',true,'componentShape','vertical-bars');
    B = cc2B(ge.cc);
    sigma_v = 0.5;
    
    nSamp = 200;
    Z = ones(nSamp,1);
    x = randn(ge.Dx,1);
        
    y_max = 0;
    for i = 1:length(g_scales)
        %G = gamrnd(2,g_scales(i),[nSamp,ge.k]); 
        if move_mean
            G = abs(normrnd(g_scales(i)*10,0.05,[nSamp,ge.k]));   
        else
            G = abs(normrnd(1,g_scales(i),[nSamp,ge.k]));
        end        
        
        g_var = var(G,0,1);
        g_mean = mean(G,1);
        g_mean_var = mean(squeeze(g_var));
        g_mean_mean = mean(squeeze(g_mean));
        
        save('bin/save_postcov.mat','G','Z');
        csm_cov = gestaltPostVCovariance(x,ge,nSamp,'leave',true,1);
        csm_corr = upperTriangleValues(corrcov(csm_cov));
        subplot(2,length(g_scales),i);
        hist(csm_corr,linspace(-1,1,100));
        title(sprintf('CSM v=%.2f m=%.2f',g_mean_var,g_mean_mean),'FontSize',16);
        yl = ylim();
        if yl(2) > y_max
            y_max = yl(2);
        end
        
        save('bin/save_msm_postcov.mat','G','Z');
        msm_cov = msmPosteriorCovariance(x,nSamp,'leave',true,ge.A,B,ge.obsVar,sigma_v,ge.g_shape,ge.g_scale,ge.z_shape,ge.z_scale,1);
        msm_corr = upperTriangleValues(corrcov(msm_cov));
        subplot(2,length(g_scales),length(g_scales) + i);
        hist(msm_corr,linspace(-1,1,100));
        title(sprintf('MSM v=%.2f m=%.2f',g_mean_var,g_mean_mean),'FontSize',16);
        yl = ylim();
        if yl(2) > y_max
            y_max = yl(2);
        end
    end
    
    for sp = 1:2*length(g_scales)
        subplot(2,length(g_scales),sp);
        xlim([-1 1]);
        ylim([0 y_max]);
    end