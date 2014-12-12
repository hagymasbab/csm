function [vsamp,gsamp,zsamp] = gestaltScheduling(stimuli,timings,models,nTrials)
    % models need to have the same dimensions and only differ in their
    % parametrisation
    if ~iscell(models)
        models = {models};
    end
    nMod = size(models,2);
    Dv = models{1}.Dv;
    k = models{1}.k;
    totalSamples = sum(timings);
    ends = cumsum(timings);
    starts = [1 ends(1:end-1)+1];
    nStim = size(stimuli,2);
    vsamp = zeros(nMod,nTrials,totalSamples,Dv);
    gsamp = zeros(nMod,nTrials,totalSamples,k);
    zsamp = zeros(nMod,nTrials,totalSamples);
    for m = 1:nMod
        fprintf('Model %d/%d ',nMod,m);
        if strcmp(models{m}.prior,'gamma') 
            g_sampler = 'gibbs-slice';
        else
            g_sampler = 'slice';
        end
        for t = 1:nTrials            
            printCounter(t,'stringVal','Trial','maxVal',nTrials,'newLine',true);
            % set initial conditions
            initZ = 0.1;
            initG = 0.1 * ones(models{m}.k,1);
            for s = 1:nStim
                % set data
                models{m}.X(1,1,:) = stimuli{s};
                %viewImage(stimuli{s});pause
                % call sampler
                [cs,~,zs] = gestaltGibbs(models{m},1,timings(s),'verbose',0,'initZ',initZ,'initG',initG,'gSampler',g_sampler);
                % store results
                vsamp(m,t,starts(s):ends(s),:) = cs(:,models{m}.k+1:end);
                gsamp(m,t,starts(s):ends(s),:) = cs(:,1:models{m}.k);
                zsamp(m,t,starts(s):ends(s)) = zs;
                % set endpoint as next initial
                initG = cs(end,1:models{m}.k)';
                initZ = zs(end);
            end
            save('gestalt_samples.mat','vsamp','gsamp','zsamp');
        end
        %fprintf('\n');
    end
end
    
    