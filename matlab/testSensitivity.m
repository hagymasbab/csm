function testSensitivity(ge,cc,nTrials,nSamples,burnin,loadSamples,plotStuff)
    % create a bunch of stimuli
    stimuli = {};
    thetas = linspace(0,120,4);
    for t = 1:length(thetas)
        lambda = 4;
        phase = 0;
        act_gr = grating(lambda,thetas(t),phase,sqrt(ge.Dx));
        stimuli{end+1} = act_gr(:);
    end
    if plotStuff
        viewImageSet(stimuli);
    end
    timings = ones(1,length(stimuli)) * (nSamples + burnin);
    if loadSamples
        load('testsens.mat');
    else
        % run scheduling
        ge.cc = cc;
        [~,gsamp,~] = gestaltScheduling(stimuli,timings,{ge},nTrials,0,true,'gibbs',false);
        save('testsens.mat','gsamp');
    end
    
    if plotStuff
        % plot g activations
        gdata = squeeze(gsamp(1,:,:,:)); % nTrials x allSamp x k
        gdata = permute(gdata,[3,1,2]); % k x nTrials x allSamp
        gdata = reshape(gdata,[1 ge.k nTrials size(gdata,3)]);
        compnums={};other={};for i=1:ge.k;compnums{end+1}=sprintf('%d',i);other{end+1}='';end; 
        plotGridSeries(gdata,cumsum(timings),other,compnums,'','C');
    end
    
end