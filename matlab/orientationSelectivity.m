function vmax = orientationSelectivity(nTrials,loadSamples,randseed,cc)

    close all;
    setrandseed(randseed);
    
    nSamples = 80;
    burnin = 20;
    
    if isempty(cc)
        Dx = 64;    
        k = 2;
        cell_idx = 41; 
        px = 0.5;
        py = 0.5;
        central_orient = 0;
    else
        Dx = size(cc{1});
        k = length(cc);
        % TODO calculate orientation, location and index of a cell
        load(sprintf('filtermatching_%d.mat',Dx));
        cell_idx = 1;
        central_orient = orients(1,1);
        px = maxX
    end
    

    %contrasts = [0.05 0.2 0.8];
    contrasts = [0.04 1];
    %contrasts = [0.5 100];
    rms_contrasts = zeros(size(contrasts));
    orient_shift = 45 * pi/180;    
    orients = central_orient;     
    for i=1:2
        orients = [central_orient-i*orient_shift orients central_orient+i*orient_shift];
    end        
        
    imSize = sqrt(Dx);
    stimuli = cell(1,length(contrasts)*length(orients));
    for o=1:length(orients)
        % create Gabor      
        [~,~,act_gabor] = gabor('theta',orients(o),'lambda',4,'width',imSize,'height',imSize,'px',px,'py',py,'Sigma',1.8);
        for z=1:length(contrasts)            
            % put Gabor on gray background
            act_stimulus = zeros(imSize) + contrasts(z) * act_gabor;
            if o==1
                rms_contrasts(z) = std(act_stimulus(:));
            end
            stimuli{1,length(orients)*(z-1)+o} = act_stimulus(:);
        end
    end
    
    viewImageSet(stimuli,'max',false);
    
    % set timings
    timings = (nSamples+burnin) * ones(1,length(contrasts)*length(orients));                
    
    if loadSamples
        load('bin/save_orient_select_samples.mat');
    else
        % create model
        if isempty(cc)
            ge = gestaltCreate('temp','Dx',Dx,'k',k,'B',1,'N',1,'filters','gabor_4or','obsVar',1, ...
                'nullComponent',false,'generateComponents',true,'generateData',false);
        else
            ge = gestaltCreate('temp','Dx',Dx,'k',k,'B',1,'N',1,'filters','OF','obsVar',1,'cc',cc, ...
                'nullComponent',false,'generateComponents',false,'generateData',false);
        end
        % run scheduling
        [vsamp,gsamp,zsamp] = gestaltScheduling(stimuli,timings,{ge},nTrials,ge.obsVar,true,'gibbs',false);
        save('bin/save_orient_select_samples.mat','vsamp','gsamp','zsamp');
    end
    
    % calculate rates
    vdata = squeeze(vsamp(1,:,:,1,cell_idx)); % ntrials x all the samples
    zdata = squeeze(zsamp(1,:,:));
    gdata = squeeze(gsamp(1,:,:,:)); % ntrials x totalsamples x k
%     figure;
%     subplot(2,1,1);
%     plot(vdata')
%     subplot(2,1,2);
%     plot(squeeze(vsamp(1,:,:,1,1))');
    %vdata = reshape(vdata,[nTrials length(contrasts) length(orients) nSamples+burnin]);
    v_split = zeros(nTrials,length(contrasts),length(orients),nSamples);
    z_split = zeros(nTrials,length(contrasts),length(orients),nSamples);
    g_split = zeros(nTrials,length(contrasts),length(orients),nSamples,k);
    for z=1:length(contrasts)
        for o=1:length(orients)
            stim_idx = (z-1)*length(orients)+o;
            start_idx = (stim_idx - 1) * (nSamples+burnin) + burnin + 1;
            end_idx = start_idx + nSamples - 1;
            v_split(:,z,o,:) = vdata(:,start_idx:end_idx);
            z_split(:,z,o,:) = zdata(:,start_idx:end_idx);
%             size(gdata(:,start_idx:end_idx,:))
%             size(g_split(:,z,o,:,:))
            g_split(:,z,o,:,:) = gdata(:,start_idx:end_idx,:);
        end
    end        
            
    %vrate = squeeze(mean(vdata(:,:,:,burnin+1:end),4));
    vrate = squeeze(mean(v_split,4)); % firing rate for each trial and stimulus nTrials x nContrast x nOrient
    vmax = squeeze(max(v_split,[],4)); % firing rate for each trial and stimulus nTrials x nContrast x nOrient
    zrate = squeeze(mean(z_split,4)); % average contrast for each trial and stimulus nTrials x nContrast x nOrient
    grate = reshape(mean(g_split,4),nTrials,length(contrasts),length(orients),k); % average g for each trial and stimulus nTrials x nContrast x nOrient x k
    size(grate)
    vvar = squeeze(var(v_split,0,4)); % variance of firing rate for each trial and stimulus
    
    figure
    scatter(zrate(:),vvar(:));
    xlabel('Estimated contrast','FontSize',16);
    ylabel('Membrane potential variance','FontSize',16);
    set(gca,'XTick',[],'YTick',[]);
    
    true_z = [];
    est_z = [];
    
    for t=1:nTrials
        for z=1:length(contrasts)
            for o=1:length(orients)
                true_z = [true_z; rms_contrasts(z)];
                est_z = [est_z; zrate(t,z,o)];
            end
        end
    end
%     rep_z = repmat(rms_contrasts,length(orients),1)';
%     true_z = repmat(rep_z(:),nTrials,1);
%     est_z = zrate(:);
    figure
    scatter(true_z,est_z);
    xlabel('RMS contrast','FontSize',16);
    ylabel('Estimated contrast','FontSize',16);
%     set(gca,'XTick',[],'YTick',[]);
    
%     v_t2t = squeeze(var(vrate,0,1)); % nContrast x nOrient
    v_t2t = [];
    rmsc = [];
    trial_avg_z = [];
    for z=1:length(contrasts)
        for o=1:length(orients)
            %v_t2t = [v_t2t; var(vmax(:,z,o))];
            v_t2t = [v_t2t; var(vrate(:,z,o))];
            rmsc = [rmsc; rms_contrasts(z)];
            trial_avg_z = [trial_avg_z; mean(zrate(:,z,o))];
        end
    end
    
%     rmsc = repmat(rms_contrasts',length(orients),1);
    figure
    scatter(rmsc,v_t2t(:));
    xlabel('RMS contrast of stimulus','FontSize',16);
    ylabel('Trial-to-trial variance of membrane potential','FontSize',16);
    %set(gca,'XTick',[],'YTick',[]);
    
%     trial_avg_z = squeeze(mean(zrate,1));
    figure
    scatter(trial_avg_z(:),v_t2t(:));
    xlabel('Estimated contrast of stimulus','FontSize',16);
    ylabel('Trial-to-trial variance of membrane potential','FontSize',16);
    
    
%     zdata = reshape(zdata,[nTrials length(contrasts) length(orients) nSamples+burnin]);
%     zrate = squeeze(mean(zdata(:,:,:,burnin+1:end),4));
    
    figure;
    % reproduce plot from s&s
    for z=1:length(contrasts)
        means = mean(vrate(:,z,:),1);
        stds = std(vrate(:,z,:),0,1);
        h = errorbar(means,stds,'LineWidth',2);
        hold on;
    end
    title('Preferred orientation of cell: 0^o','FontSize',16);
    xlabel('Stimulus orientation','FontSize',16)
    ylabel(sprintf('Firing rate, mean and std of %d trials',nTrials),'FontSize',16)
    ticlabels = cell(1,length(orients));
    for o=1:length(orients)
        ticlabels{o} = sprintf('%d',orients(o)*180/pi);
    end
    set(gca,'XTick',1:5,'XTickLabels',ticlabels,'FontSize',16,'Ytick',[]);
    
    figure;
    for z=1:length(contrasts)
        means = mean(zrate(:,z,:),1);
        stds = std(zrate(:,z,:),0,1);
        h = errorbar(means,stds,'LineWidth',2);        
        hold on;
    end
    title('Preferred orientation of cell: 0^o','FontSize',16);
    xlabel('Stimulus orientation','FontSize',16)
    ylabel(sprintf('Estimated contrast, mean and std of %d trials',nTrials),'FontSize',16)    
    %set(gca,'XTick',1:5,'XTickLabels',ticlabels,'FontSize',16,'Ytick',[]);
    
    gcontr = [];
    gstds = [];
    for z=1:length(contrasts)
        actg_m = [];
        actg_s = [];
        for o=1:length(orients)
            size(grate(:,z,o,:))
            g = reshape(grate(:,z,o,:),nTrials,k); % ntrials x k
            %size(g)
            actg_m = [actg_m; mean(g,1)]; % o x k
            actg_s = [actg_s; std(g,0,1)]; % o x k
        end
        %size(actg_m)
        gcontr = [gcontr; mean(actg_m,1)]; % z x k
        gstds = [gstds; mean(actg_s,1)]; % z x k
    end
    
    figure;
    for kk=1:k
        means = gcontr(:,kk);
        stds = gstds(:,kk);
        h = errorbar(means,stds,'LineWidth',2);        
        hold on;
    end
    title('G','FontSize',16);
end