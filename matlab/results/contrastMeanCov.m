function contrastMeanCov(Dv,randseed,loadStuff,plotStuff)

    setrandseed(randseed);   
    imdim = sqrt(Dv);
    % load filter set
    load(sprintf('filters_gabor_4or_%d.mat',Dv));
    
    % define covariance components
    cc{1} = eye(Dv);
    if Dv == 256
        filter_ids = [12 16 30];
    elseif Dv == 64
        filter_ids = [6 8 15];
    end
    for f = 1:length(filter_ids)
        for ff = f+1:length(filter_ids)
            i = filter_ids(f);
            j = filter_ids(ff);
            cc{1}(i,j) = 1;
            cc{1}(j,i) = 1;
            cc{1}(i,i) = cc{1}(i,i) + 1;
            cc{1}(j,j) = cc{1}(j,j) + 1;
        end        
    end    
    
    %filter_ids = [2*imdim+2*4 2*imdim+3*4 2*imdim+4*4];
    %filter_ids = [146-4*16+8 146 146+4*16-8];
    if Dv == 256
        filter_ids = [146-2*16+4 146 146+2*16-4];
    elseif Dv == 64
        filter_ids = [36-2*8+4 36 36+2*8-4];
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
    %B = randn(Dv,2);  
    B = zeros(Dv,2);
    for f = 1:length(filter_ids)
        B(filter_ids(f),1) = 1;
    end
    sigma_v = 0.5;
    
    ge = gestaltCreate('temp','Dx',Dv,'k',length(cc),'filters','gabor_4or','obsVar',0.5,'cc',cc, ...
        'g_shape',2,'g_scale',2,'z_shape',2,'z_scale',2,'N',1,'generateComponents',false,'generateData',false);

    % construct base stimulus
    %x_base = zeros(Dv,1);
%     x_base = A * randn(Dv,1);
%     for f = 1:length(filter_ids)
%         x_base = x_base + A(:,filter_ids(f));
%     end
    x_base = gestaltAncestralSample(ge,[1;0],1)';
    %viewImage(x_base);
    
    % create contrast-adjusted stimuli
    %contrasts = [0.01 0.1 0.5 1 2 5 10 20 100];
    contrasts = [0.5 2];
    
    x_rms = zeros(length(contrasts),1);
%     corr_c = zeros(length(contrasts),1);
%     corr_c_ort_priv = zeros(length(contrasts),1);
%     corr_c_ort = zeros(length(contrasts),1);
%     corr_c_priv = zeros(length(contrasts),1);
%     corr_m = zeros(length(contrasts),1);
%     corr_m_priv = zeros(length(contrasts),1);
%     corr_m_ort = zeros(length(contrasts),1);
%     corr_m_ort_priv = zeros(length(contrasts),1);
%     
%     var_c_gest = zeros(length(contrasts),1);
%     var_c_priv = zeros(length(contrasts),1);
%     var_m_gest = zeros(length(contrasts),1);
%     var_m_priv = zeros(length(contrasts),1);
%     
%     cell1 = filter_ids(1);
%     cell2 = filter_ids(2);
%     cell3 = filter_ids(3);
       
    nSamp = 50;
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
            [covc,gsamp,zsamp] = posteriorCovariances(x_act,ge,nSamp,randseed,false);
            [covm,gsampm,zsampm] = msmPosteriorCovariance(x_act,nSamp,randseed,false,ge.A,B,ge.obsVar,sigma_v,ge.g_shape,ge.g_scale,ge.z_shape,ge.z_scale);
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
    
        
%         % get posterior covariances for each stimuli    
%         [covc,covm] = posteriorCovariances(x_act,A,0.5,0.5,2,2,2,2,cc,B,200,randseed);
%         % get correlations and variances
%         corrc = corrcov(covc);
%         corrm = corrcov(covm);
%         
%         corr_c(c) = corrc(cell1,cell2);
%         corr_c_ort(c) = corrc(cell1,cell3);
%         corr_c_priv(c) = corrc(12,16);
%         corr_c_ort_priv(c) = corrc(12,30);        
%         
%         corr_m(c) = corrm(cell1,cell2);    
%         corr_m_ort(c) = corrm(cell1,cell3);
%         corr_m_priv(c) = corrm(12,16);
%         corr_m_ort_priv(c) = corrm(12,30);
%         
%         var_c_gest(c) = covc(cell1,cell1);
%         var_m_gest(c) = covm(cell1,cell1);
%         var_c_priv(c) = covc(10,10);
%         var_m_priv(c) = covm(10,10);
%     end
    % plot results
    
%     subplot(2,1,1)
%     plot(log(x_rms),abs([corr_c corr_c_ort corr_c_priv corr_c_ort_priv corr_m corr_m_ort corr_m_priv corr_m_ort_priv]),'LineWidth',2);
%     legend({'C ovl cov','C ort cov','C ovl priv','C ort priv','M ovl cov','M ort cov','M ovl priv','M ort priv'})
%     xlim([0 max(log(x_rms))]);    
%     ylabel('Correlation magnitude','FontSize',16);
%     subplot(2,1,2)
%     plot(log(x_rms),[var_c_gest var_c_priv var_m_gest var_m_priv],'LineWidth',2);
%     legend({'C cov','C priv','M cov','M priv'});
%     xlim([0 max(log(x_rms))]);    
%     xlabel('log-RMS contrast of stimulus','FontSize',16);
%     ylabel('Variance','FontSize',16);
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