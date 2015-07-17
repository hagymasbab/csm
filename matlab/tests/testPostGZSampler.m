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
        numcols = 3;
        
        mom0_mean_g = zeros(nSN,ge.k);
        mom0_std_g = zeros(nSN,ge.k);
        mom0_mean_z = zeros(nSN,1);
        mom0_std_z = zeros(nSN,1);
        
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
            mom0_trials_g = zeros(nTrials,ge.k);
            mom0_trials_z = zeros(nTrials,1);
            for t=1:nTrials
                for kk = 1:ge.k
                    mom0_trials_g(t,kk) = sampleMode(squeeze(gsamples{i}(t,:,kk)),100);
                end
                mom0_trials_z(t,1) = sampleMode(squeeze(zsamples{i}(t,:,1)),100);
            end
            
            mom1_trials_g = squeeze(mean(gsamples{i},2)); % nTrials x k
            mom1_trials_z = squeeze(mean(zsamples{i},2)); % nTrials x 1
            mom2_trials_g = squeeze(var(gsamples{i},0,2)); % nTrials x k
            mom2_trials_z = squeeze(var(zsamples{i},0,2)); % nTrials x 1
            
            % calculate mean and std for these moments over trials
            mom0_mean_g(i,:) = mean(mom0_trials_g,1);
            mom0_std_g(i,:) = std(mom0_trials_g,0,1);    
            mom0_mean_z(i,:) = mean(mom0_trials_z,1);
            mom0_std_z(i,:) = std(mom0_trials_z,0,1);  
            
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
        
        subplot(2,numcols,1)
        barwitherr(mom0_std_g,mom0_mean_g);
        set(gca,'XTickLabel',snlabels,'FontSize',16);
        xlabel('# of samples');
        ylabel({'Mean and STD';'zeroth sample moment of G'});
        title(sprintf('Dv=%d, k=%d, nTrials=%d, burn-in=%d',ge.Dv,ge.k,nTrials,burnin));
        subplot(2,numcols,numcols+1)
        barwitherr(mom0_std_z,mom0_mean_z);
        set(gca,'XTickLabel',snlabels,'FontSize',16);
        xlabel('# of samples');
        ylabel({'Mean and STD';'zeroth sample moment of Z'});
        
        subplot(2,numcols,2)
        barwitherr(mom1_std_g,mom1_mean_g);
        set(gca,'XTickLabel',snlabels,'FontSize',16);
        xlabel('# of samples');
        ylabel({'Mean and STD';'first sample moment of G'});        
        subplot(2,numcols,numcols+2)
        barwitherr(mom1_std_z,mom1_mean_z);
        set(gca,'XTickLabel',snlabels,'FontSize',16);
        xlabel('# of samples');
        ylabel({'Mean and STD';'first sample moment of Z'});
        
        subplot(2,numcols,3)    
        barwitherr(mom2_std_g,mom2_mean_g);
        set(gca,'XTickLabel',snlabels,'FontSize',16);
        xlabel('# of samples');
        ylabel({'Mean and STD';'second sample moment of G'});
        subplot(2,numcols,numcols+3)
        barwitherr(mom2_std_z,mom2_mean_z);
        set(gca,'XTickLabel',snlabels,'FontSize',16);
        xlabel('# of samples');
        ylabel({'Mean and STD';'second sample moment of Z'});
        
        subplot(2,numcols,1);
        title(sprintf('Dv=%d, k=%d, nTrials=%d, burn-in=%d',ge.Dv,ge.k,nTrials,burnin));
        subplot(2,numcols,2);
        gstr_parts = [];
        for kk=1:ge.k
            gstr_parts = [gstr_parts sprintf(',%.2f',ge.G(1,kk))];
        end
        gstr_parts = gstr_parts(2:end);
        title(['True G=[' gstr_parts sprintf('] Z=%.2f',ge.Z(1,1))]);
    end
end
