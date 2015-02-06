function orientationSelectivity(nTrials,loadSamples,randseed)

    close all;
    setrandseed(randseed);
    Dx = 64;
    imSize = sqrt(Dx);
    k = 2;
    nSamples = 80;
    burnin = 20;

    %contrasts = [0.05 0.2 0.8];
    contrasts = [0.5 2 3 5];
    orient_shift = 45 * pi/180;
    central_orient = 0;
    orients = [central_orient-2*orient_shift central_orient-orient_shift central_orient central_orient+orient_shift central_orient+2*orient_shift];        
    
    px = 0.5;
    py = 0.5;
    
    % calculate cell id    
    cell_idx = 41; 
    
    stimuli = cell(1,length(contrasts)*length(orients));
    for o=1:length(orients)
        % create Gabor      
        [~,~,act_gabor] = gabor('theta',orients(o),'lambda',4,'width',imSize,'height',imSize,'px',px,'py',py,'Sigma',1.8);
        for z=1:length(contrasts)
            % put Gabor on gray background
            act_stimulus = zeros(imSize) + contrasts(z) * act_gabor;
            stimuli{1,length(orients)*(z-1)+o} = act_stimulus(:);
        end
    end
    
    viewImageSet(stimuli,'max',false);
    
    % set timings
    timings = (nSamples+burnin) * ones(1,length(contrasts)*length(orients));        
    
    % create model
    ge = gestaltCreate('temp','Dx',Dx,'k',k,'B',1,'N',1,'filters','gabor_4or','obsVar',1,'z_shape',1, ...
        'nullComponent',false,'generateComponents',true,'generateData',false);
    
    if loadSamples
        load('orient_select_samples.mat');
    else
        % run scheduling
        [vsamp,gsamp,zsamp] = gestaltScheduling(stimuli,timings,{ge},nTrials,ge.obsVar,true,'gibbs');
        save('orient_select_samples.mat','vsamp','gsamp','zsamp');
    end
    
    % calculate rates
    vdata = squeeze(vsamp(1,:,:,1,cell_idx)); % ntrials x all the samples
    zdata = squeeze(zsamp(1,:,:));
%     figure;
%     subplot(2,1,1);
%     plot(vdata')
%     subplot(2,1,2);
%     plot(squeeze(vsamp(1,:,:,1,1))');
    %vdata = reshape(vdata,[nTrials length(contrasts) length(orients) nSamples+burnin]);
    v_split = zeros(nTrials,length(contrasts),length(orients),nSamples);
    z_split = zeros(nTrials,length(contrasts),length(orients),nSamples);
    for z=1:length(contrasts)
        for o=1:length(orients)
            stim_idx = (z-1)*length(orients)+o;
            start_idx = (stim_idx - 1) * (nSamples+burnin) + burnin + 1;
            end_idx = start_idx + nSamples - 1;
            v_split(:,z,o,:) = vdata(:,start_idx:end_idx);
            z_split(:,z,o,:) = zdata(:,start_idx:end_idx);
        end
    end        
            
    %vrate = squeeze(mean(vdata(:,:,:,burnin+1:end),4));
    vrate = squeeze(mean(v_split,4));
    zrate = squeeze(mean(z_split,4));
    vvar = squeeze(var(v_split,0,4));
    figure
    scatter(zrate(:),vvar(:));
    xlabel('Estimated contrast','FontSize',16);
    ylabel('Membrane potential variance','FontSize',16);
    set(gca,'XTick',[],'YTick',[]);
    
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
    set(gca,'XTick',1:5,'XTickLabels',ticlabels,'FontSize',16,'Ytick',[]);
    
end