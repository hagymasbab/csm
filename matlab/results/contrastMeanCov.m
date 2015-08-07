function contrastMeanCov(Dv,randseed,loadStuff,plotStuff,target_acceptance,nSamp,contrasts,stimulus)

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
    if strcmp(stimulus,'msm')
        %x_base = msmGenerate(1,'leave',ge.A,B,ge.obsVar,sigma_v,ge.g_shape,ge.g_scale,ge.z_shape,ge.z_scale)';
        x_base = msmAncestralSample(genG,genZ,'leave',ge.A,B,ge.obsVar,sigma_v)';
    elseif strcmp(stimulus,'csm')
        x_base = gestaltAncestralSample(ge,genG,genZ)';    
    elseif strcmp(stimulus,'natural')
        load(sprintf('patches_vanhateren_%d.mat',Dv));
        x_base = patchDB(:,randi(size(patchDB,2),1));
        base_rms = std(x_base);
        contrasts = contrasts * (1/base_rms);
    else
        error('Ivalid stimulus generation choice: %s',stimulus);
    end
    figure;viewImage(x_base);    
    
    x_rms = zeros(length(contrasts),1);
    nCont = length(contrasts);
    filterCatPerc = 0.1;
        
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
            
            cat_c = categorizeFilters(x_act,A,ge.cc,'CSM',filterCatPerc,false);
            cat_m = categorizeFilters(x_act,A,B,'MSM',filterCatPerc,false);
            cat_g = categorizeFilters(x_act,A,C_gsm,'GSM',filterCatPerc,false);
%             cat_c = 0;
%             cat_m = 0;
%             cat_g = 0;
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
        figure;
        plotCatCorr(contrasts,corrmats,csm_cats,'CSM',false,3,1);
        plotCatCorr(contrasts,msm_corrmats,msm_cats,'MSM',false,3,2);
        plotCatCorr(contrasts,gsm_corrmats,gsm_cats,'GSM',true,3,3);
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

function plotCatCorr(xrms,corrmats,cats,modelName,x_axis,nModels,number)
    %figure;
    nCont = length(xrms);
    nCat = length(cats{1}.categories);
%     nRow = floor(sqrt(nCat));
%     nCol = ceil(nCat/nRow);

    cat_corr_std = ones(nCont,nCat);
    cat_corr_mean = ones(nCont,nCat);
    labels = {};
    legend_names = cats{1}.categories;
    
    for cont = 1:nCont
        labels{end+1} = sprintf('%.2f',xrms(cont));
        for cat = 1:nCat
    %         subplot(nRow,nCol,cat);   
            act_assign = cats{cont}.category_assignments{cat};
            nPairs = size(act_assign,1);
            if cont == 1
                %legend_names{cat} = [legend_names{cat} sprintf(', N=%d',nPairs)];
                legend_names{cat} = [legend_names{cat} sprintf('%d',nPairs)];
            end
            corrvalues = zeros(nPairs,1);
            for pair = 1:nPairs
                corrvalues(pair) = abs(corrmats{cont}(act_assign(pair,1),act_assign(pair,2)));
            end
            cat_corr_std(cont,cat) = std(corrvalues,0,1) / nPairs;
            cat_corr_mean(cont,cat) = mean(corrvalues,1);
        end
    end
    
%     barwitherr(cat_corr_std,cat_corr_mean)    
%     if x_axis
%         xlabel('Z_{RMS}','FontSize',16);
%         set(gca,'XTickLabel',labels,'FontSize',16);
%     else
%         set(gca,'XTickLabel',{});
%     end
%     legend(legend_names,'FontSize',16);
    
    if nCat > 20
        %subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0.001], [0 0.025], [0 0.01]);
        fLen = floor(nCat / 2);
        subplot(nModels*2,1,(number-1)*2+1);
        barwitherr(cat_corr_std(:,1:fLen)',cat_corr_mean(:,1:fLen)');
        set(gca,'XTickLabel',legend_names(1:fLen),'FontSize',16);
        %legend(labels,'FontSize',16);
        
        %lh=findall(gcf,'tag','legend');
        %set(lh,'location','northeastoutside');
        title(sprintf('%s',modelName),'FontSize',16);
        
        subplot(nModels*2,1,(number-1)*2+2);
        barwitherr(cat_corr_std(:,fLen+1:end)',cat_corr_mean(:,fLen+1:end)');
        set(gca,'XTickLabel',legend_names(fLen+1:end),'FontSize',16);
        %legend(labels,'FontSize',16);                
    else
        subplot(nModels,1,number);
        barwitherr(cat_corr_std',cat_corr_mean')    
        if x_axis
            %xlabel('Z_{RMS}','FontSize',16);
            set(gca,'XTickLabel',legend_names,'FontSize',16);
        else
            set(gca,'XTickLabel',{});
        end
        legend(labels,'FontSize',16);
    
    
        ylabel({'V post. corr. magn.';'mean, SEM'},'FontSize',16);
        lh=findall(gcf,'tag','legend');
        set(lh,'location','northeastoutside');
        title(sprintf('Model = %s',modelName),'FontSize',16);
    end
end
