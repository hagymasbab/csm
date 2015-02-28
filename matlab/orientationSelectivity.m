function vmax = orientationSelectivity(nTrials,loadSamples,randseed,cc,plotStuff)

    close all;
    setrandseed(randseed);
    
    nSamples = 150;
    burnin = 100;
%     nSamples = 20;
%     burnin = 5;
    
    if isempty(cc)
        Dx = 64;    
        k = 2;
        cell_idx = 41; 
        px = 0.5;
        py = 0.5;
        central_orient = 0;
    else
        Dx = size(cc{1},1);
        k = length(cc);
        % calculate orientation, location and index of a cell
        load(sprintf('filtermatching_%d.mat',Dx));
        %load('filtermatching_576.mat');
        cell_idx = 3;
        central_orient = orients(cell_idx,1);
        %if central_orient > 90; central_orient = 180 - central_orient;end;
        px = maxX(1,cell_idx) / sqrt(Dx);
        py = maxY(1,cell_idx) / sqrt(Dx);
    end
    

    %contrasts = [0.05 0.2 0.8];
    contrasts = [0.05 10 50];
    %contrasts = [0.5 100];
    rms_contrasts = zeros(size(contrasts));
    stepnum = 2;
    orient_shift = (45/stepnum);    
    stim_orients = central_orient;     
    for i=1:stepnum
        stim_orients = [central_orient-i*orient_shift stim_orients central_orient+i*orient_shift];
    end        
    %stim_orients = [central_orient-(stepnum+1)*orient_shift stim_orients];
        
    imSize = sqrt(Dx);
    stimuli = cell(1,length(contrasts)*length(stim_orients));
    for o=1:length(stim_orients)
        % create Gabor      
        [~,~,act_gabor] = gabor('theta',circ_ang2rad(stim_orients(o)),'lambda',4,'width',imSize,'height',imSize,'px',px,'py',py,'Sigma',1.8);
        for z=1:length(contrasts)            
            % put Gabor on gray background
            act_stimulus = zeros(imSize) + contrasts(z) * act_gabor;
            if o==1
                rms_contrasts(z) = std(act_stimulus(:));
            end
            stimuli{1,length(stim_orients)*(z-1)+o} = act_stimulus(:);
        end
    end        
    
    %viewImageSet(stimuli,'max',false);
    %pause
    
    % set timings
    timings = (nSamples+burnin) * ones(1,length(contrasts)*length(stim_orients));                
    
    % create model
    if isempty(cc)
        ge = gestaltCreate('temp','Dx',Dx,'k',k,'B',1,'N',1,'filters','gabor_4or','obsVar',1, ...
            'nullComponent',false,'generateComponents',true,'generateData',false);
    else
        ge = gestaltCreate('temp','Dx',Dx,'k',k,'B',1,'N',1,'filters','OF','obsVar',1,'cc',cc, ...
            'nullComponent',false,'generateComponents',false,'generateData',false);
    end
    %figure;viewImage(ge.A(:,cell_idx)','useMax',true);
    
    if loadSamples
        load('bin/save_orient_select_samples.mat');
    else        
        % run scheduling
        [vsamp,gsamp,zsamp] = gestaltScheduling(stimuli,timings,{ge},nTrials,ge.obsVar,true,'gibbs',false);
        save('bin/save_orient_select_samples.mat','vsamp','gsamp','zsamp');
    end
 
    if plotStuff
        vdata = reshape(vsamp(1,:,:,1,cell_idx),[nTrials sum(timings)]); 
        zdata = reshape(zsamp(1,:,:),[nTrials sum(timings)]);
        gdata = reshape(gsamp(1,:,:,:),[nTrials sum(timings) k]);

        v_split = zeros(nTrials,length(contrasts),length(stim_orients),nSamples);
        z_split = zeros(nTrials,length(contrasts),length(stim_orients),nSamples);
        g_split = zeros(nTrials,length(contrasts),length(stim_orients),nSamples,k);
        for z=1:length(contrasts)
            for o=1:length(stim_orients)
                stim_idx = (z-1)*length(stim_orients)+o;
                start_idx = (stim_idx - 1) * (nSamples+burnin) + burnin + 1;
                end_idx = start_idx + nSamples - 1;
                v_split(:,z,o,:) = vdata(:,start_idx:end_idx);
                z_split(:,z,o,:) = zdata(:,start_idx:end_idx);
                g_split(:,z,o,:,:) = gdata(:,start_idx:end_idx,:);
            end
        end        

        vrate = reshape(mean(v_split,4),[nTrials length(contrasts) length(stim_orients)]); % firing rate for each trial and stimulus 
        vvar = reshape(var(v_split,0,4),[nTrials length(contrasts) length(stim_orients)]); % variance of firing rate for each trial and stimulus
        vmax = reshape(max(v_split,[],4),[nTrials length(contrasts) length(stim_orients)]); % firing rate for each trial and stimulus
        zrate = reshape(mean(z_split,4),[nTrials length(contrasts) length(stim_orients)]); % average contrast for each trial and stimulus

        figure
        scatter(zrate(:),vvar(:));
        xlabel('Estimated contrast','FontSize',16);
        ylabel('Membrane potential variance','FontSize',16);
        set(gca,'XTick',[],'YTick',[]);

        true_z = [];
        est_z = [];

        for t=1:nTrials
            for z=1:length(contrasts)
                for o=1:length(stim_orients)
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
    
        figure
        scatter(true_z,vvar(:));
        xlabel('RMS contrast','FontSize',16);
        ylabel('Membrane potential variance','FontSize',16);
        set(gca,'XTick',[],'YTick',[]);
        
        zbar_mean = [];
        zbar_std = [];
        xticklabels = {};
        vars = vvar(:);        
        for z=1:length(rms_contrasts)
            actvvars = vars(true_z==rms_contrasts(z));
            zbar_mean = [zbar_mean; mean(actvvars)];
            zbar_std = [zbar_std; std(actvvars)];
            xticklabels{end+1} = sprintf('%.2f',rms_contrasts(z));
        end
        figure
        barwitherr(zbar_std,zbar_mean);
        set(gca,'XTickLabel',xticklabels);
        xlabel('RMS contrast','FontSize',16);
        ylabel('Membrane potential variance','FontSize',16);

    % %     v_t2t = squeeze(var(vrate,0,1)); % nContrast x nOrient
    %     v_t2t = [];
    %     rmsc = [];
    %     trial_avg_z = [];
    %     for z=1:length(contrasts)
    %         for o=1:length(orients)
    %             %v_t2t = [v_t2t; var(vmax(:,z,o))];
    %             v_t2t = [v_t2t; var(vrate(:,z,o))];
    %             rmsc = [rmsc; rms_contrasts(z)];
    %             trial_avg_z = [trial_avg_z; mean(zrate(:,z,o))];
    %         end
    %     end
    %     
    %     figure
    %     scatter(rmsc,v_t2t(:));
    %     xlabel('RMS contrast of stimulus','FontSize',16);
    %     ylabel('Trial-to-trial variance of membrane potential','FontSize',16);
    %     %set(gca,'XTick',[],'YTick',[]);
    %     
    %     figure
    %     scatter(trial_avg_z(:),v_t2t(:));
    %     xlabel('Estimated contrast of stimulus','FontSize',16);
    %     ylabel('Trial-to-trial variance of membrane potential','FontSize',16);


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
        title(sprintf('Matched orientation of cell: %.2f^o',central_orient),'FontSize',16);
        xlabel('Stimulus orientation','FontSize',16)
        ylabel(sprintf('Firing rate, mean and std of %d trials',nTrials),'FontSize',16)
        ticlabels = cell(1,length(stim_orients));
        for o=1:length(stim_orients)
            ticlabels{o} = sprintf('%0.2f',stim_orients(o));
        end
        set(gca,'XTick',1:length(stim_orients),'XTickLabel',ticlabels,'FontSize',16,'Ytick',[]);

        figure;
        for z=1:length(contrasts)
            means = mean(zrate(:,z,:),1);
            stds = std(zrate(:,z,:),0,1);
            h = errorbar(means,stds,'LineWidth',2);        
            hold on;
        end
        %title('Preferred orientation of cell: 0^o','FontSize',16);
        xlabel('Stimulus orientation','FontSize',16)
        ylabel(sprintf('Estimated contrast, mean and std of %d trials',nTrials),'FontSize',16)    
        %set(gca,'XTick',1:5,'XTickLabels',ticlabels,'FontSize',16,'Ytick',[]);

    %     grate = reshape(mean(g_split,4),nTrials,length(contrasts),length(stim_orients),k); % average g for each trial and stimulus nTrials x nContrast x nOrient x k
    %     gcontr = [];
    %     gstds = [];
    %     for z=1:length(contrasts)
    %         actg_m = [];
    %         actg_s = [];
    %         for o=1:length(orients)
    %             %size(grate(:,z,o,:))
    %             g = reshape(grate(:,z,o,:),nTrials,k); % ntrials x k
    %             %size(g)
    %             actg_m = [actg_m; mean(g,1)]; % o x k
    %             actg_s = [actg_s; std(g,0,1)]; % o x k
    %         end
    %         %size(actg_m)
    %         gcontr = [gcontr; mean(actg_m,1)]; % z x k
    %         gstds = [gstds; mean(actg_s,1)]; % z x k
    %     end
    %     
    %     figure;
    %     for kk=1:k
    %         means = gcontr(:,kk);
    %         stds = gstds(:,kk);
    %         h = errorbar(means,stds,'LineWidth',2);        
    %         hold on;
    %     end
    %     title('G','FontSize',16);
    end
end