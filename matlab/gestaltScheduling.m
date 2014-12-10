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
    nStim = size(stimuli,2);
    vsamp = zeros(nMod,nTrials,totalSamples,Dv);
    gsamp = zeros(nMod,nTrials,totalSamples,k);
    zsamp = zeros(nMod,nTrials,totalSamples);
    for m = 1:nMod
        for t = 1:nTrials
            
end
    
    