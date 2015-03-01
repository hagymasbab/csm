function testSensitivity(ge,cc,nTrials,allSamples,burnin,loadSamples,plotStuff)
    % create a bunch of stimuli
    stimuli = {};
    %thetas = linspace(0,120,4);
%     thetas = [90 140];
%     for t = 1:length(thetas)
%         lambda = 10;
%         phase = 0.5;
%         act_gr = grating(lambda,thetas(t),phase,sqrt(ge.Dx));
%         %act_gr = act_gr + grating(lambda,thetas(t)+90,phase,sqrt(ge.Dx));
%         stimuli{end+1} = act_gr(:);
%     end
    load img_stone
    stimuli{1} = img(:);
    load img_seeds
    stimuli{2} = img(:);
%     if plotStuff
%         viewImageSet(stimuli);
%         pause
%     end
    
    timings = ones(1,length(stimuli)) * allSamples;
    if loadSamples
        load('bin/save_testsens.mat');
    else
        % run scheduling
        ge.cc = cc;
        [vsamp,gsamp,zsamp] = gestaltScheduling(stimuli,timings,{ge},nTrials,0,true,'gibbs',false);
        save('bin/save_testsens.mat','gsamp','vsamp','zsamp','stimuli');
    end
    
    if plotStuff
        close all;
        % plot g activations
        gdata = squeeze(gsamp(1,:,:,:)); % nTrials x allSamp x k
        gdata = permute(gdata,[3,1,2]); % k x nTrials x allSamp
        gdata = reshape(gdata,[1 ge.k nTrials size(gdata,3)]);
        compnums={};other={};for i=1:ge.k;compnums{end+1}=sprintf('%d',i);other{end+1}='';end; 
        %plotGridSeries(gdata,cumsum(timings),other,compnums,'','C');
        
        xticpos = cumsum(timings) - timings(1)/2;
        responses = zeros(length(stimuli),ge.k);
        resp_stds = zeros(length(stimuli),ge.k);
        for i=1:ge.k
            xticlab = {};
            for j=1:length(stimuli)
                start_idx = (j-1)*timings(1) + burnin + 1;
                end_idx = j*timings(1);
                act_g = reshape(gdata(1,i,:,start_idx:end_idx),nTrials,allSamples-burnin); % nTrials x nSamples
                trial_means = mean(act_g,2);
                mean_response = mean(trial_means);
                resp_stds(j,i) = std(trial_means);
                responses(j,i) = mean_response;
                xticlab{end+1} = sprintf('mu = %.2f',mean_response);
            end
            %subplot(ge.k,1,i);
            %set(gca,'XTick',xticpos,'XTickLabel',xticlab,'Ytick',[]);
        end
        
        for j=1:length(stimuli)
            figure;
            viewImage(stimuli{j},'useMax',true);
        end
        
        figure();
        load cmp_graybars;
        ymax = 0;
        for j=1:length(stimuli)
            subplot(length(stimuli),1,j);
            barwitherr(resp_stds(j,:)'/sqrt(nTrials),responses(j,:)');
            xlim([0 ge.k+1]);
            set(gca,'XTickLabel',compnums,'YTick',[],'FontSize',16);
            if j==1
                title('Activation of components');
            elseif j==length(stimuli)
                xlabel('Component #','FontSize',16);
            end
            colormap(cmp_graybars);
            ylims = ylim();
            if ylims(2) > ymax
                ymax = ylims(2);
            end
            %title(sprintf('Stimulus %d',j));
            %subplot(length(stimuli),2,j*2-1);
            %viewImage(stimuli{j},'useMax',true);
        end
        
        for j=1:length(stimuli)
            subplot(length(stimuli),1,j);
            ylim([0 ymax]);
        end
    end
    
end