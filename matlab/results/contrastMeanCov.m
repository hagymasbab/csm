function contrastMeanCov(Dv,randseed,loadStuff,plotStuff,target_acceptance,nSamp,contrasts,genFromMSM)

    setrandseed(randseed);               
    imdim = sqrt(Dv);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%% DEFINE MODEL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % load filter set
    load(sprintf('filters_gabor_4or_%d.mat',Dv));
    
    % define covariance components
%     cc{1} = eye(Dv);
%     if Dv == 256
%         filter1_ids = [12 16 30];
%     elseif Dv == 64
%         filter1_ids = [6 8 15];
%     elseif Dv == 576
%         filter1_ids = [24 32 60];
%     end
%     for f = 1:length(filter1_ids)
%         for ff = f+1:length(filter1_ids)
%             i = filter1_ids(f);
%             j = filter1_ids(ff);
%             cc{1}(i,j) = 1;
%             cc{1}(j,i) = 1;
%             cc{1}(i,i) = cc{1}(i,i) + 1;
%             cc{1}(j,j) = cc{1}(j,j) + 1;
%         end        
%     end    
%     
%     if Dv == 256
%         filter_ids = [146-2*16+4 146 146+2*16-4];
%     elseif Dv == 64
%         filter_ids = [36-2*8+4 36 36+2*8-4];
%     elseif Dv == 576
%         filter_ids = [292-2*24+4 292 292+2*24-4];
%     end
%     cc{2} = eye(Dv);
%     for f = 1:length(filter_ids)
%         for ff = f+1:length(filter_ids)
%             i = filter_ids(f);
%             j = filter_ids(ff);
%             cc{2}(i,j) = 1;
%             cc{2}(j,i) = 1;
%             cc{2}(i,i) = cc{2}(i,i) + 1;
%             cc{2}(j,j) = cc{2}(j,j) + 1;
%         end        
%     end                
%     
    tempGe.Dv = Dv;
    tempGe.Dx = Dv;
    cc = gestaltCovariances(tempGe,2,'nullComponent',false,'method','vertical-bars');
    
%     viewImageSet(cc);
%     figure
    
    % define mean components
%     B = zeros(Dv,2);
%     for f = 1:length(filter1_ids)
%         B(filter1_ids(f),1) = 1;
%     end
%     for f = 1:length(filter_ids)
%         B(filter_ids(f),2) = 1;
%     end
    B = cc2B(cc);
    
    sigma_v = 0.5;
    
    % define residual component for CSM
    cc{3} = sigma_v * eye(Dv);
    
    % define GSM prior covariance
%     randM = randn(Dv);
%     C_gsm = randM * randM';
%     C_gsm = C_gsm / (max(C_gsm(:)));
    C_gsm = componentSum(1,cc(1:2));
%     C_gsm = sigma_v * eye(Dv);
    
    ge = gestaltCreate('temp','Dx',Dv,'k',length(cc),'filters','gabor_4or','obsVar',0.7,'cc',cc, ...
        'g_shape',1,'g_scale',0.1,'z_shape',2,'z_scale',2,'N',1,'generateComponents',false,'generateData',false,'nullComponent',length(cc)==3);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% DEFINE STIMULUS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    genG = zeros(ge.k,1);
    genG(2) = 5;
    genZ = 1;
    if genFromMSM
        %x_base = msmGenerate(1,'leave',ge.A,B,ge.obsVar,sigma_v,ge.g_shape,ge.g_scale,ge.z_shape,ge.z_scale)';
        x_base = msmAncestralSample(genG,genZ,'leave',ge.A,B,ge.obsVar,sigma_v)';
    else
        x_base = gestaltAncestralSample(ge,genG,genZ)';    
    end
    %viewImage(x_base);    
    
    x_rms = zeros(length(contrasts),1);
    nCont = length(contrasts);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% CALCULATE POSTERIORS %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if loadStuff
        load('save_contmeancov.mat');
    else        
        corrmats = {};
        covmats = {};
        gsamps = {};
        zsamps = {};
        csm_cats = {};
        
        msm_corrmats = {};
        msm_covmats = {};        
        msm_gsamps = {};
        msm_zsamps = {};        
        msm_cats = {};
        
        gsm_corrmats = {};
        gsm_covmats = {};        
        gsm_zmoments = {};        
        gsm_cats = {};
        
        for c = 1:nCont
            x_act = contrasts(c) * x_base;
            
%             cat_c = categorizeFilters(x_act,A,ge.cc,'CSM');
%             cat_m = categorizeFilters(x_act,A,B,'MSM');
%             cat_g = categorizeFilters(x_act,A,C_gsm,'GSM');
            cat_c = 0;
            cat_m = 0;
            cat_g = 0;
            csm_cats{end+1} = cat_c;
            msm_cats{end+1} = cat_m;
            gsm_cats{end+1} = cat_g;
            
            [covc,gsamp,zsamp] = gestaltPostVCovariance(x_act,ge,nSamp,randseed,false,target_acceptance);
            [covm,gsampm,zsampm] = msmPosteriorCovariance(x_act,nSamp,randseed,false,ge.A,B,ge.obsVar,sigma_v,ge.g_shape,ge.g_scale,ge.z_shape,ge.z_scale,target_acceptance);
            [~,covg,zmeang,zstdg] = gsmPosteriorV(x_act,ge.A,C_gsm,sqrt(ge.obsVar),ge.z_shape,ge.z_scale,30);
            
            covmats{end+1} = covc;
            msm_covmats{end+1} = covm;
            gsm_covmats{end+1} = covg;
            corrmats{end+1} = corrcov(covc);
            msm_corrmats{end+1} = corrcov(covm);
            gsm_corrmats{end+1} = corrcov(covg);            
            
            gsamps{end+1} = gsamp;
            zsamps{end+1} = zsamp;
            msm_gsamps{end+1} = gsampm;
            msm_zsamps{end+1} = zsampm;
            gsm_zmoments{end+1} = [zmeang zstdg];
        end
        save('bin/save_contmeancov.mat','covmats','corrmats','csm_cats','gsamps','zsamps', ...
                                        'msm_covmats','msm_corrmats','msm_cats','msm_gsamps','msm_zsamps', ...
                                        'gsm_covmats','gsm_corrmats','gsm_cats','gsm_zmoments');
    end           
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if plotStuff    
        close all;
        plotResults(contrasts,x_base,covmats,corrmats,gsamps,zsamps,'CSM');
        plotResults(contrasts,x_base,msm_covmats,msm_corrmats,msm_gsamps,msm_zsamps,'MSM');
        plotResults(contrasts,x_base,gsm_covmats,gsm_corrmats,[],gsm_zmoments,'GSM');
        %plotCatCorr(contrasts,corrmats,csm_cats,'CSM');
        %plotCatCorr(contrasts,msm_corrmats,msm_cats,'MSM');
    end
    
end

function plotResults(contrasts,x_base,covmats,corrmats,gsamps,zsamps,modelName)
    figure;
    nCont = length(contrasts);
    nCol = max(nCont,3);
    nRow = 3;
    x_rms = zeros(nCont,1);
    
    if ~strcmp(modelName,'GSM')
        k = size(gsamps{1},2);
    else
        k = 0;
    end
    g_means = zeros(nCont,k);
    g_stds = zeros(nCont,k);
    z_means = zeros(nCont,1);
    z_stds = zeros(nCont,1);
    vvar_means = zeros(nCont,1);
    vvar_stds = zeros(nCont,1);
    g_labels = {};
    hist_y_max = 0;
    for c = 1:nCont
        x_act = contrasts(c) * x_base;        
        x_rms(c) = std(x_act(:));
        
        subplot(nRow,nCol,c);
        viewImage(covmats{c});
        xlabel(sprintf('Covariance, Z_{RMS} = %.2f',x_rms(c)),'FontSize',16);
        if c==1
            title(sprintf('Model = %s',modelName),'FontSize',16);
        end
        
        subplot(nRow,nCol,nCol+c);
        hist(upperTriangleValues(corrmats{c}),linspace(-1,1,100));
        yl = ylim();
        if yl(2) > hist_y_max
            hist_y_max = yl(2);
        end
        %viewImage(corrmats{c});
        xlabel(sprintf('Correlations, Z_{RMS} = %.2f',x_rms(c)),'FontSize',16);
        
        if ~strcmp(modelName,'GSM')
            g_means(c,:) = mean(gsamps{c},1);
            g_stds(c,:) = std(gsamps{c},0,1);
            z_means(c,1) = mean(zsamps{c},1);
            z_stds(c,1) = std(zsamps{c},0,1);
        else
            z_means(c,1) = zsamps{c}(1);
            z_stds(c,1) = zsamps{c}(2);
        end
        
        g_labels{end+1} = sprintf('%.2f',x_rms(c));
        vvar_means(c,1) = mean(diag(covmats{c}));
        vvar_stds(c,1) = std(diag(covmats{c}));
    end
    
    for c = 1:nCont
        subplot(nRow,nCol,nCol+c);
        ylim([0 hist_y_max]);
        xlim([-1 1]);
    end
    
    if ~strcmp(modelName,'GSM')
        subplot(nRow,nCol,2*nCol+1);
        barwitherr(g_stds,g_means);
        xlabel('Z_{RMS}','FontSize',16);
        ylabel('G posterior mean and std','FontSize',16);
        set(gca,'XTickLabel',g_labels,'FontSize',16);
        ylim([0 max(g_means(:))+max(g_stds(:))+0.1])
    end
    
    subplot(nRow,nCol,2*nCol+2);
    barwitherr(z_stds,z_means);
    xlabel('Z_{RMS}','FontSize',16);
    ylabel('Z posterior mean and std','FontSize',16);
    set(gca,'XTickLabel',g_labels,'FontSize',16);
    
    subplot(nRow,nCol,2*nCol+3);
    barwitherr(vvar_stds,vvar_means);
    xlabel('Z_{RMS}','FontSize',16);
    ylabel('V post. variance mean and std','FontSize',16);
    set(gca,'XTickLabel',g_labels,'FontSize',16);
    ylim([0 max(vvar_means)+max(vvar_stds)+0.1])
end

function plotCatCorr(contrasts,corrmats,cats,modelName)
    figure;
    nCont = length(contrasts);
    for c = 1:nCont
        
    end
end
