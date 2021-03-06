function plotICSampleVariances(filename,z)
%function plotICSampleVariances(allsamp,within_var,trial_var,within_cov,central_field,A,k,prestim,poststim,nullComp)
    
    close all;
    load(filename);
    
    if nullComp
        model_k = k+1;
    else
        model_k = k;
    end        
    
    orstrings = {'|','/','-','\\'};
    within_var_cells = squeeze(var(allsamp(:,z,:,prestimSamp+1:end-poststimSamp,model_k+central_field'),0,3));
    %within_var_cells = squeeze(within_var(:,z,:,model_k+central_field'));
    within_var_means = squeeze(mean(within_var_cells,2));
    within_var_stds = squeeze(std(within_var_cells,0,2));
    %test
    %within_var_means = within_var_means(1:3,:);
    %within_var_stds = within_var_stds(1:3,:);
    
    trial_var_cells = squeeze(trial_to_trial_variance(:,z,:,model_k+central_field'));
    trial_var_means = squeeze(mean(trial_var_cells,2));
    trial_var_stds = squeeze(std(trial_var_cells,0,2));               
    
    figure();
    subplot(2,1,1)
    barwitherr(within_var_stds,within_var_means);
    title('Within-trial variance of samples averaged over trials','FontSize',16);
    legend(orstrings,'FontSize',16);
    subplot(2,1,2)
    barwitherr(trial_var_stds,trial_var_means);
    title('Trial-to-trial variance of samples averaged over sampling steps','FontSize',16);
    
    within_cov_cells = squeeze(within_trial_covariance(:,z,:,model_k+central_field,model_k+1:end));
    within_cov_means = squeeze(mean(within_cov_cells,2)); % -> nStim x nOrient x Dv
    nStim = size(allsamp,1);
    nOrient = size(central_field,1);
    Dv = size(allsamp,5) - model_k;
    cov_images = cell(nStim,nOrient);   
    covtitles = {};
    for stim = 1:nStim
        for ori = 1:nOrient
            coeffs = reshape(within_cov_means(stim,ori,:),Dv,1);
            cov_images{stim,ori} = A * coeffs;
            covtitles{end+1} = sprintf('Stim: %s RF: %s',orstrings{stim},orstrings{ori});
        end
    end
    
    figure();
    viewImageSet(cov_images,'max',false,'titles',covtitles);
        
    vsamples = squeeze(allsamp(:,z,:,:,model_k+central_field')); % nStim x nTrial x nSamp x nOrient
    figure();
    for stim = 1:nStim
        for ori = 1:nOrient
            subplot(nStim,nOrient,(stim-1)*nOrient+ori);
            actdata = squeeze(vsamples(stim,:,:,ori)); % nTrial x nSamp
            plot(actdata');
            ylim([-1 2]);
            xlim([1,size(allsamp,4)]);
            hold on;
            plot(mean(actdata)','LineWidth',3);
            title(sprintf('Stim: %s RF: %s',orstrings{stim},orstrings{ori}),'FontSize',16);
        end
    end
    
    gsamples = squeeze(allsamp(:,z,:,:,1:model_k)); % nStim x nTrial x nSamp x k
    gtitles = orstrings;
    if nullComp
        gtitles{end+1} = '0';
    end
    figure();
    for stim = 1:nStim
        for gest = 1:model_k
            subplot(nStim,model_k,(stim-1)*model_k+gest);
            actdata = squeeze(gsamples(stim,:,:,gest)); % nTrial x nSamp
            plot(actdata');
            ylim([0 1]);
            xlim([1,size(allsamp,4)]);
            hold on;
            plot(mean(actdata)','LineWidth',3);
            title(sprintf('Stim: %s Gestalt: %s',orstrings{stim},gtitles{gest}),'FontSize',16);
        end
    end
    
    zsamples = squeeze(zsamp(:,z,:,:)); % nStim x nTrial x nSamp
    figure();
    for stim = 1:nStim
        subplot(nStim,1,stim);
        actdata = squeeze(zsamples(stim,:,:)); % nTrial x nSamp
        plot(actdata');
        ylim([0 5]);
        xlim([1,size(allsamp,4)]);
        hold on;
        plot(mean(actdata)','LineWidth',3);
        title(sprintf('Z for stim: %s',orstrings{stim}),'FontSize',16);
    end
    
    actvsamples = squeeze(allsamp(:,z,:,:,model_k+activated_units')); % nStim x nTrial x nSamp x nOrient
    figure();
    for stim = 1:nStim
        for unit = 1:length(activated_units)
            subplot(nStim,length(activated_units),(stim-1)*length(activated_units)+unit);
            actdata = squeeze(actvsamples(stim,:,:,unit)); % nTrial x nSamp
            plot(actdata');
            ylim([-1 2]);
            xlim([1,size(allsamp,4)]);
            hold on;
            plot(mean(actdata)','LineWidth',3);
            title(sprintf('Stim: %s',orstrings{stim}),'FontSize',16);
        end
    end
end