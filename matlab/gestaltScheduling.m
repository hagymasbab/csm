function [vsamp,gsamp,zsamp] = gestaltScheduling(stimuli,timings,models,nTrials,obsNoise,reset,sampler,breakdownSamples)
    if sum(timings~=timings(1)) > 0 && breakdownSamples
        error('cannot break down samples if timings are not equal');
    end
    % models need to have the same dimensions and only differ in their
    % parametrisation
    if ~iscell(models)
        models = {models};
    end
    if ~iscell(stimuli)
        stimuli = mat2cell(stimuli,ones(1,size(stimuli,1)));
    end
    nMod = size(models,2);
    Dv = models{1}.Dv;
    B = models{1}.B;
    ks = zeros(1,nMod);
    for m = 1:nMod
        ks(m) = models{m}.k;
    end
    k = max(ks);
    totalSamples = sum(timings);
    nStim = size(stimuli,2);
    if ~breakdownSamples
        ends = cumsum(timings);
        starts = [1 ends(1:end-1)+1];        
        vsamp = zeros(nMod,nTrials,totalSamples,B,Dv);
        gsamp = zeros(nMod,nTrials,totalSamples,k);
        zsamp = zeros(nMod,nTrials,totalSamples);
    else
        vsamp = zeros(nMod,nTrials,nStim,timings(1),B,Dv);
        gsamp = zeros(nMod,nTrials,nStim,timings(1),k);
        zsamp = zeros(nMod,nTrials,nStim,timings(1));
    end
    for m = 1:nMod
        fprintf('Model %d/%d ',nMod,m);
        if strcmp(models{m}.prior,'gamma') 
            g_sampler = 'gibbs-slice';
            %g_sampler = 'slice';
        else
            g_sampler = 'slice';
        end
        for t = 1:nTrials            
            printCounter(t,'stringVal','Trial','maxVal',nTrials,'newLine',true);
            % set initial conditions
            initZ = 1;
            initG = 0.5 * ones(models{m}.k,1);
            for s = 1:nStim
                % set data
                actstim = zeros(models{m}.B,models{m}.Dx);
                % TODO repmat
                for bb = 1 : models{m}.B
                    actstim(bb,:) = stimuli{s}(:) + obsNoise * randn(size(stimuli{s}));
                end
                models{m}.X(1,:,:) = actstim;
                %viewImage(models{m}.X(1,1,:));pause
                
                % call sampler
                if strcmp(sampler,'gibbs')
                    [vs,gs,zs,~] = gestaltGibbs(models{m},1,timings(s),'verbose',0,'initZ',initZ,'initG',initG,'gSampler',g_sampler);
                elseif strcmp(sampler,'hamilton')
                    [vs,gs,zs,~] = gestaltHamiltonian(models{m},actstim,timings(s),'sampler','nuts','burnin',timings(s));
                end
                
                % store results
                if ~breakdownSamples
                    vsamp(m,t,starts(s):ends(s),:,:) = vs;
                    gsamp(m,t,starts(s):ends(s),1:models{m}.k) = gs;
                    zsamp(m,t,starts(s):ends(s)) = zs;
                else
                    vsamp(m,t,s,:,:,:) = vs;
                    gsamp(m,t,s,:,1:models{m}.k) = gs;
                    zsamp(m,t,s,:) = zs;
                end
                
                % set endpoint as next initial
                if ~reset                    
                    initG = gs(end,:)';
                    initZ = zs(end);
                end
            end
            save('gestalt_samples.mat','vsamp','gsamp','zsamp');
        end
        %fprintf('\n');
    end
end
    
    