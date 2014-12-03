function plotICSampleVariances(allsamp,within_var,trial_var,within_cov,central_field,A,k)
    
    close all;
    model_k = k+1;
    
    z = 2;
    orstrings = {'|','/','-','\\'};

    within_var_cells = squeeze(within_var(:,z,:,model_k+central_field'));
    within_var_means = squeeze(mean(within_var_cells,2));
    within_var_stds = squeeze(std(within_var_cells,0,2));
    %test
    %within_var_means = within_var_means(1:3,:);
    %within_var_stds = within_var_stds(1:3,:);
    
    trial_var_cells = squeeze(trial_var(:,z,:,model_k+central_field'));
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
    
    within_cov_cells = squeeze(within_cov(:,z,:,model_k+central_field,model_k+1:end));
    within_cov_means = squeeze(mean(within_cov_cells,2)); % -> nStim x nOrient x Dv
    nStim = size(allsamp,1);
    nOrient = size(central_field,1);
    Dv = size(allsamp,5) - model_k;
    cov_images = cell(nStim,nOrient);    
    for stim = 1:nStim
        for ori = 1:nOrient
            coeffs = reshape(within_cov_means(stim,ori,:),Dv,1);
            cov_images{stim,ori} = A * coeffs;
        end
    end
    
    figure();
    viewImageSet(cov_images);
        
    vsamples = squeeze(allsamp(:,z,:,:,model_k+central_field')); % nStim x nTrial x nSamp x nOrient
    figure();
    for stim = 1:nStim
        for ori = 1:nOrient
            subplot(nStim,nOrient,(stim-1)*nOrient+ori);
            actdata = squeeze(vsamples(stim,:,:,ori)); % nTrial x nSamp
            plot(actdata');
            ylim([-5 10]);
            xlim([1,size(allsamp,4)]);
            hold on;
            plot(mean(actdata)','LineWidth',3);
            title(sprintf('Stim: %s RF: %s',orstrings{stim},orstrings{ori}),'FontSize',16);
        end
    end
    
    gsamples = squeeze(allsamp(:,z,:,:,1:model_k)); % nStim x nTrial x nSamp x k
    gtitles = orstrings;
    gtitles{end+1} = '0';
    figure();
    for stim = 1:nStim
        for gest = 1:model_k
            subplot(nStim,model_k,(stim-1)*model_k+gest);
            actdata = squeeze(gsamples(stim,:,:,gest)); % nTrial x nSamp
            plot(actdata');
            ylim([0 12]);
            xlim([1,size(allsamp,4)]);
            hold on;
            plot(mean(actdata)','LineWidth',3);
            title(sprintf('Stim: %s Gestalt: %s',orstrings{stim},gtitles{gest}),'FontSize',16);
        end
    end
end