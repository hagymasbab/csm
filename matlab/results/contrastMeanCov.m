function contrastMeanCov(Dv,randseed,loadStuff,plotStuff,target_acceptance,nSamp,contrasts,genFromMSM)

    setrandseed(randseed);   
    imdim = sqrt(Dv);
    % load filter set
    load(sprintf('filters_gabor_4or_%d.mat',Dv));
    
    % define covariance components
    cc{1} = eye(Dv);
    if Dv == 256
        filter1_ids = [12 16 30];
    elseif Dv == 64
        filter1_ids = [6 8 15];
    elseif Dv == 576
        filter1_ids = [24 32 60];
    end
    for f = 1:length(filter1_ids)
        for ff = f+1:length(filter1_ids)
            i = filter1_ids(f);
            j = filter1_ids(ff);
            cc{1}(i,j) = 1;
            cc{1}(j,i) = 1;
            cc{1}(i,i) = cc{1}(i,i) + 1;
            cc{1}(j,j) = cc{1}(j,j) + 1;
        end        
    end    
    
    if Dv == 256
        filter_ids = [146-2*16+4 146 146+2*16-4];
    elseif Dv == 64
        filter_ids = [36-2*8+4 36 36+2*8-4];
    elseif Dv == 576
        filter_ids = [292-2*24+4 292 292+2*24-4];
    end
    cc{2} = eye(Dv);
    for f = 1:length(filter_ids)
        for ff = f+1:length(filter_ids)
            i = filter_ids(f);
            j = filter_ids(ff);
            cc{2}(i,j) = 1;
            cc{2}(j,i) = 1;
            cc{2}(i,i) = cc{2}(i,i) + 1;
            cc{2}(j,j) = cc{2}(j,j) + 1;
        end        
    end    
%     viewImageSet(cc);
%     figure
    
    % define mean components
    B = zeros(Dv,2);
    for f = 1:length(filter1_ids)
        B(filter1_ids(f),1) = 1;
    end
    for f = 1:length(filter_ids)
        B(filter_ids(f),2) = 1;
    end
    sigma_v = 0.5;
    
    ge = gestaltCreate('temp','Dx',Dv,'k',length(cc),'filters','gabor_4or','obsVar',0.1,'cc',cc, ...
        'g_shape',1,'g_scale',0.1,'z_shape',2,'z_scale',2,'N',1,'generateComponents',false,'generateData',false);

    genG = [0;5];
    genZ = 1;
    if genFromMSM
        %x_base = msmGenerate(1,'leave',ge.A,B,ge.obsVar,sigma_v,ge.g_shape,ge.g_scale,ge.z_shape,ge.z_scale)';
        x_base = msmAncestralSample(genG,genZ,'leave',ge.A,B,ge.obsVar,sigma_v)';
    else
        x_base = gestaltAncestralSample(ge,genG,genZ)';    
    end
    viewImage(x_base);    
    
    x_rms = zeros(length(contrasts),1);
    nCont = length(contrasts);
    
    if loadStuff
        load('save_contmeancov.mat');
    else        
        corrmats = {};
        gsamps = {};
        zsamps = {};
        msm_corrmats = {};
        msm_gsamps = {};
        msm_zsamps = {};
        for c = 1:nCont
            x_act = contrasts(c) * x_base;            
            [covc,gsamp,zsamp] = posteriorCovariances(x_act,ge,nSamp,randseed,false,target_acceptance);
            [covm,gsampm,zsampm] = msmPosteriorCovariance(x_act,nSamp,randseed,false,ge.A,B,ge.obsVar,sigma_v,ge.g_shape,ge.g_scale,ge.z_shape,ge.z_scale,target_acceptance);
            corrmats{end+1} = corrcov(covc);
            msm_corrmats{end+1} = corrcov(covm);
            gsamps{end+1} = gsamp;
            zsamps{end+1} = zsamp;
            msm_gsamps{end+1} = gsampm;
            msm_zsamps{end+1} = zsampm;
        end
        save('bin/save_contmeancov.mat','corrmats','gsamps','zsamps','msm_corrmats','msm_gsamps','msm_zsamps');
    end           
    
    if plotStuff    
        close all;
        plotResults(contrasts,x_base,corrmats,gsamps,zsamps,'CSM');
        plotResults(contrasts,x_base,msm_corrmats,msm_gsamps,msm_zsamps,'MSM');
    end
    
end

function plotResults(contrasts,x_base,corrmats,gsamps,zsamps,modelName)
    figure;
    nCont = length(contrasts);
    x_rms = zeros(nCont,1);
    
    k = size(gsamps{1},2);
    g_means = zeros(nCont,k);
    g_stds = zeros(nCont,k);
    z_means = zeros(nCont,1);
    z_stds = zeros(nCont,1);
    g_labels = {};
    for c = 1:nCont
        x_act = contrasts(c) * x_base;        
        x_rms(c) = std(x_act(:));
        subplot(2,nCont,c);
        viewImage(corrmats{c});
        xlabel(sprintf('%s post. corr. of v, Z_{RMS} = %.2f',modelName,x_rms(c)),'FontSize',16);
        g_means(c,:) = mean(gsamps{c},1);
        g_stds(c,:) = std(gsamps{c},0,1);
        z_means(c,1) = mean(zsamps{c},1);
        z_stds(c,1) = std(zsamps{c},0,1);
        g_labels{end+1} = sprintf('%.2f',x_rms(c));
    end
    subplot(2,nCont,nCont+1);
    barwitherr(real(g_stds),real(g_means));
    xlabel('Z_{RMS}','FontSize',16);
    ylabel('G posterior mean and std','FontSize',16);
    set(gca,'XTickLabel',g_labels,'FontSize',16);
    subplot(2,nCont,nCont+2);
    barwitherr(real(z_stds),real(z_means));
    xlabel('Z_{RMS}','FontSize',16);
    ylabel('Z posterior mean and std','FontSize',16);
    set(gca,'XTickLabel',g_labels,'FontSize',16);
end