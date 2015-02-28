function responseCorrelations(samplefile,cc,burnin)
    close all;

    % load samples from the file
    load(sprintf('save_testsens_%d.mat',samplefile));
    % vsamp: 1 x nTrial x nStim*nSamp x 1 x Dv
    nStim = length(stimuli);
    nSamp = size(vsamp,3) / nStim;
    nTrials = size(vsamp,2);
    Dx = size(vsamp,5);
    k = size(gsamp,4);    
    
    vdata = reshape(vsamp(1,:,:,1,:),[nTrials nStim*nSamp Dx]);
    gdata = reshape(gsamp(1,:,:,:),[nTrials nStim*nSamp k]);
    zdata = reshape(zsamp(1,:,:),[nTrials nStim*nSamp]);
    v_split = zeros(nTrials,nStim,nSamp-burnin,Dx);
    g_split = zeros(nTrials,nStim,nSamp-burnin,k);
    z_split = zeros(nTrials,nStim,nSamp-burnin);
    for s=1:nStim
        v_split(:,s,:,:) = vdata(:,(s-1)*nSamp+1+burnin:s*nSamp,:);
        g_split(:,s,:,:) = gdata(:,(s-1)*nSamp+1+burnin:s*nSamp,:);
        z_split(:,s,:) = zdata(:,(s-1)*nSamp+1+burnin:s*nSamp);
    end
    vrate = reshape(mean(v_split,3),[nTrials nStim Dx]);
    grate = reshape(mean(g_split,3),[nTrials nStim k]);
    zrate = reshape(mean(z_split,3),[nTrials nStim]);
    vmean = reshape(mean(v_split,1),[nStim nSamp-burnin Dx]);
    grate_mean = reshape(mean(grate),nStim,k);        
    
    % calculate correlation matrix of v-samples for each trial
    wt_corrmats = cell(nStim,nTrials);
    t2t_corrmats = cell(1,nStim);
    wt_mean_corrmats = cell(1,nStim); % correlations of the mean responses
    average_corrmats = cell(1,nStim); % mean correlations of the responses in each trial
    for s=1:nStim
        for t=1:nTrials
            act_v = squeeze(v_split(t,s,:,:)); % nSamp-burnin and Dx should not be 1
            wt_corrmats{s,t} = corr(act_v);
        end
        act_rate = squeeze(vrate(:,s,:)); % nTrials and Dx should not be 1
        t2t_corrmats{1,s} = corr(act_rate);
        act_mean = squeeze(vmean(s,:,:)); % nSamp-burnin and Dx should not be 1
        wt_mean_corrmats{1,s} = corr(act_mean);
        average_corrmats{1,s} = componentSum(1/k,wt_corrmats(s,:));
    end
%     figure;
%     viewImageSet(wt_mean_corrmats);
%     figure;
%     viewImageSet(t2t_corrmats);
    
    cv = componentSum(1/k,cc);
    cv_post = cell(1,nStim);
    cr_post = cell(1,nStim);
    stimlabels = {};
    for s=1:nStim
        cv_post{s} = componentSum(grate_mean(s,:)',cc);
        cr_post{s} = corrcov(cv_post{s});
        stimlabels{s} = sprintf('Stim %d',s);
    end
    
    load('cmp_graybars.mat');
    cmp = summer;
    
    %figure;
    bar(zrate);
    xlim([0 nTrials+1]);
    xlabel('Trial #','FontSize',16);
    ylabel('Estimated contrast','FontSize',16);
    set(gca,'FontSize',16);
    legend(stimlabels,'FontSize',16);
    colormap(cmp);
    
    %[~,cell1,cell2] = maxNElements(nodiag(corrcov(cc{2})),1);
    %[~,cell1,cell2] = maxNElements(nodiag(cc{2})-nodiag(cv),1);
    %[~,cell1,cell2] = maxNElements(nodiag(corrcov(cc{2}))-nodiag(corrcov(cv)),1);
    [~,cell1,cell2] = maxNElements(nodiag(wt_mean_corrmats{2}),1);
    vars_covs_prior = zeros(k,4);
    for i = 1:k
        corrm = corrcov(cc{i});
        vars_covs_prior(i,1) = cc{i}(cell1,cell1);
        vars_covs_prior(i,2) = cc{i}(cell2,cell2);
        vars_covs_prior(i,3) = cc{i}(cell1,cell2); 
        vars_covs_prior(i,4) = corrm(cell1,cell2);        
    end
    figure;
    subplot(4,1,1);
    bar(vars_covs_prior);
    xlim([0 k+1]);
    legend({'Variance cell 1','Variance cell 2','Covariance','Correlation'},'FontSize',16)
    xlabel('Component #','FontSize',16);
    title('Prior statistics per component');
    set(gca,'FontSize',16);
    colormap(cmp);
    
    %figure;
    subplot(4,1,2);
    bar(grate_mean');
    xlim([0 k+1]);
    legend(stimlabels,'FontSize',16)
    xlabel('Component #','FontSize',16);
    ylabel('g sample mean','FontSize',16);
    title('Component activation');
    set(gca,'FontSize',16);
    colormap(cmp);
    
    vars_covs_posterior = zeros(nStim,4);
    for s=1:nStim            
        vars_covs_posterior(s,1) = cv_post{s}(cell1,cell1);
        vars_covs_posterior(s,2) = cv_post{s}(cell2,cell2);
        vars_covs_posterior(s,3) = cv_post{s}(cell1,cell2);
        vars_covs_posterior(s,4) = cr_post{s}(cell1,cell2);
    end
    %figure;
    subplot(4,1,3);
    bar(vars_covs_posterior);
    xlim([0 nStim+1]);
    legend({'Variance cell 1','Variance cell 2','Covariance','Correlation'},'FontSize',16)
    xlabel('Stimulus #','FontSize',16);
    title('Posterior statistics');
    set(gca,'FontSize',16);
    colormap(cmp);
    
    corrs = [wt_mean_corrmats{1}(cell1,cell2) wt_mean_corrmats{2}(cell1,cell2); ...
    t2t_corrmats{1}(cell1,cell2) t2t_corrmats{2}(cell1,cell2); ...
    average_corrmats{1}(cell1,cell2) average_corrmats{2}(cell1,cell2)];
    %figure;
    subplot(4,1,4);
    bar(corrs);
    xlim([0 size(corrs,1)+1]);
    set(gca,'XTickLAbel',{'Mean resp. corr.','T2T corr','Mean corr. of resp'},'FontSize',16)
    legend(stimlabels,'FontSize',16);
    title('Sample statistics');
    set(gca,'FontSize',16);
    colormap(cmp);
    
end