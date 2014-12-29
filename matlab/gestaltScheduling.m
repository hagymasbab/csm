function [vsamp,gsamp,zsamp] = gestaltScheduling(stimuli,timings,models,nTrials,obsNoise,reset)
    % models need to have the same dimensions and only differ in their
    % parametrisation
    if ~iscell(models)
        models = {models};
    end
    nMod = size(models,2);
    Dv = models{1}.Dv;
    B = models{1}.B;
    k = models{1}.k;
    totalSamples = sum(timings);
    ends = cumsum(timings);
    starts = [1 ends(1:end-1)+1];
    nStim = size(stimuli,2);
    vsamp = zeros(nMod,nTrials,totalSamples,B,Dv);
    gsamp = zeros(nMod,nTrials,totalSamples,k);
    zsamp = zeros(nMod,nTrials,totalSamples);
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
            initZ = 0.1;
            initG = 0.1 * ones(models{m}.k,1);
            for s = 1:nStim
                % set data
                actstim = zeros(models{m}.B,models{m}.Dx);
                % TODO repmat
                for bb = 1 : models{m}.B
                    actstim(bb,:) = stimuli{s} + obsNoise * randn(size(stimuli{s}));
                end
                models{m}.X(1,:,:) = actstim;
                %viewImage(models{m}.X(1,1,:));pause
                
                % call sampler
                [cs,~,zs] = gestaltGibbs(models{m},1,timings(s),'verbose',0,'initZ',initZ,'initG',initG,'gSampler',g_sampler);
                
                % store results
                actlength = ends(s) - starts(s) + 1;
                vsamp(m,t,starts(s):ends(s),:,:) = reshape(cs(:,models{m}.k+1:end),[actlength B Dv]);
                gsamp(m,t,starts(s):ends(s),:) = cs(:,1:models{m}.k);
                zsamp(m,t,starts(s):ends(s)) = zs;
                
                % set endpoint as next initial
                if ~reset                    
                    initG = cs(end,1:models{m}.k)';
                    initZ = zs(end);
                end
            end
            save('gestalt_samples.mat','vsamp','gsamp','zsamp');
        end
        %fprintf('\n');
    end
end
    
    