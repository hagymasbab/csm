function testPostGZSampler(randseed,Dv,k,sampleNums,nTrials,loadSamples,plotStuff)
    % see if the estimations of the first two moments converge as we
    % increase sample size
    
    setrandseed(randseed);
        
    ge = gestaltCreate('temp','Dx',Dv,'k',k,'filters','gabor_4or','obsVar',0.5,'N',1, ...
        'g_shape',2,'g_scale',2,'z_shape',2,'z_scale',2,'N',1,'generateComponents',true,'generateData',true);
    
    ge.Z
    ge.G
    
    x = reshape(ge.X(1,1,:),ge.Dx,1);
    burnin = 10;
    
    nSN = length(sampleNums);
    if loadSamples
        load('save_testgzsamp.mat');
    else
        gsamples = cell(1,nSN);
        zsamples = cell(1,nSN);
        for i = 1:nSN
            gtr = zeros(nTrials,sampleNums(i),ge.k);
            ztr = zeros(nTrials,sampleNums(i),1);
            for t = 1:nTrials
                [gtr(t,:,:),ztr(t,:,:)] = gestaltHamiltonianGZ(x,ge,sampleNums(i),burnin,'shuffle');
            end
            gsamples{i} = gtr;
            zsamples{i} = ztr;
        end
        save('bin/save_testgzsamp.mat','gsamples','zsamples');
    end
    
    if plotStuff
        mom1_mean_g = zeros(nSN,ge.k);
        mom1_std_g = zeros(nSN,ge.k);
        mom1_mean_z = zeros(nSN,1);
        mom1_std_z = zeros(nSN,1);
        
        mom2_mean_g = zeros(nSN,ge.k);
        mom2_std_g = zeros(nSN,ge.k);
        mom2_mean_z = zeros(nSN,1);
        mom2_std_z = zeros(nSN,1);
        
        snlabels = {};
        for i = 1:nSN
            % calculate sample moments for each trial
            mom1_trials_g = squeeze(mean(gsamples{i},2)); % nTrials x k
            mom1_trials_z = squeeze(mean(zsamples{i},2)); % nTrials x 1
            mom2_trials_g = squeeze(var(gsamples{i},0,2)); % nTrials x k
            mom2_trials_z = squeeze(var(zsamples{i},0,2)); % nTrials x 1
            
            % calculate mean and std for these moments over trials
            mom1_mean_g(i,:) = mean(mom1_trials_g,1);
            mom1_std_g(i,:) = std(mom1_trials_g,0,1);    
            mom1_mean_z(i,:) = mean(mom1_trials_z,1);
            mom1_std_z(i,:) = std(mom1_trials_z,0,1);   
            
            mom2_mean_g(i,:) = mean(mom2_trials_g,1);
            mom2_std_g(i,:) = std(mom2_trials_g,0,1);    
            mom2_mean_z(i,:) = mean(mom2_trials_z,1);
            mom2_std_z(i,:) = std(mom2_trials_z,0,1);   
            
            snlabels{end+1} = sprintf('%d',sampleNums(i));
        end
        subplot(2,2,1)        
        barwitherr(mom1_std_g,mom1_mean_g);
        set(gca,'XTickLabel',snlabels,'FontSize',16);
        xlabel('# of samples');
        ylabel('Mean and STD of first sample moment of G');
        title(sprintf('Dv=%d, k=%d, nTrials=%d, burn-in=%d',ge.Dv,ge.k,nTrials,burnin));
        subplot(2,2,2)
        barwitherr(mom1_std_z,mom1_mean_z);
        set(gca,'XTickLabel',snlabels,'FontSize',16);
        xlabel('# of samples');
        ylabel('Mean and STD of first sample moment of Z');
        
        subplot(2,2,3)        
        barwitherr(mom2_std_g,mom2_mean_g);
        set(gca,'XTickLabel',snlabels,'FontSize',16);
        xlabel('# of samples');
        ylabel('Mean and STD of second sample moment of G');        
        subplot(2,2,4)
        barwitherr(mom2_std_z,mom2_mean_z);
        set(gca,'XTickLabel',snlabels,'FontSize',16);
        xlabel('# of samples');
        ylabel('Mean and STD of second sample moment of Z');
    end
end