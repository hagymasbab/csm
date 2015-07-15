function contrastMeanCov(randseed,loadStuff,plotStuff)

    setrandseed(randseed);
    Dv = 256;
    imdim = sqrt(Dv);
    % load filter set
    load(sprintf('filters_gabor_4or_%d.mat',Dv));
    
    % define covariance components
    cc{1} = eye(Dv);
    filter_ids = [12 16 30];
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
    filter_ids = [146-2*16+4 146 146+2*16-4];
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
%     B = zeros(Dv,2);
%     for f = 1:length(filter_ids)
%         B(filter_ids(f),1) = 1;
%     end
    
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
    contrasts = [0.05 0.1 0.5 1 5 10];
    
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
        for c = 1:nCont
            x_act = contrasts(c) * x_base;            
            [covc,gsamp,zsamp] = posteriorCovariances(x_act,ge,nSamp,randseed,false);            
            corrmats{end+1} = corrcov(covc);
            gsamps{end+1} = gsamp;
            zsamps{end+1} = zsamp;
        end
        save('bin/save_contmeancov.mat','corrmats','gsamps','zsamps');
    end           
    
    if plotStuff        
        g_means = zeros(nCont,ge.k);
        g_stds = zeros(nCont,ge.k);
        z_means = zeros(nCont,1);
        z_stds = zeros(nCont,1);
        g_labels = {};
        for c = 1:nCont
            x_act = contrasts(c) * x_base;        
            x_rms(c) = std(x_act(:));
            subplot(2,nCont,c);
            viewImage(corrmats{c});
            xlabel(sprintf('Posterior correlation of v, Z_{RMS} = %.2f',x_rms(c)),'FontSize',16);
            g_means(c,:) = mean(gsamps{c},1);
            g_stds(c,:) = std(gsamps{c},0,1);
            z_means(c,1) = mean(zsamps{c},1);
            z_stds(c,1) = std(zsamps{c},0,1);
            g_labels{end+1} = sprintf('%.2f',x_rms(c));
        end
        subplot(2,nCont,nCont+1);
        barwitherr(g_stds,g_means);
        xlabel('Z_{RMS}','FontSize',16);
        ylabel('G posterior mean and std','FontSize',16);
        set(gca,'XTickLabel',g_labels,'FontSize',16);
        subplot(2,nCont,nCont+2);
        barwitherr(z_stds,z_means);
        xlabel('Z_{RMS}','FontSize',16);
        ylabel('Z posterior mean and std','FontSize',16);
         set(gca,'XTickLabel',g_labels,'FontSize',16);
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