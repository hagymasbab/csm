function gestaltPriorCorrelations(nTrials,timings,appendTo,calculation)
    close all;
    Dx = 1024;    
    B = 10;
    cumulative_timings = cumsum(timings);
    
    nOrient = 4; % this is true for gabor_4or_32.mat only
    shift = nOrient/2; 
    rfsinarow = sqrt(Dx)/shift;    
        
    z_shape = 1;
    z_scale = 0.1;
    sigma = 0.1;
    g_mean = 0.1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CREATE MODELS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ic_idx = floor(Dx / (nOrient*2) + rfsinarow/2);
    sh1 = 3; % single shift 
    sh2 = 5; % double shift
    locations = [ic_idx-(sh2*rfsinarow)+sh2 ic_idx-(sh1*rfsinarow)+sh1 ic_idx+(sh1*rfsinarow)-sh1 ic_idx+(sh2*rfsinarow)-sh2];
    
    % create covariance components
    filterlists = [ ([ic_idx locations] - 1) .* nOrient + 2; ([ic_idx locations] - 10) .* nOrient + 3 ];
    cc = filterList2Components(filterlists,true,Dx);
    k = size(cc,2); 

    models = {};
    model_names = {};
    
    g_shape = 1;
    ge1 = gestaltCreate('temp','Dx',Dx,'k',k,'B',B,'N',1,'nullComponent',true,'filters','gabor_4or_32.mat','cc',cc, ...
        'obsVar',sigma,'z_shape',z_shape,'z_scale',z_scale, ...
        'prior','gamma','g_shape',g_shape,'g_scale',g_mean/g_shape,'null_shape',g_shape,'null_scale',g_mean/g_shape);
    models{end+1} = ge1;
    model_names{end+1} = 'gamma1';
    
    g_shape = 2;
    ge2 = gestaltCreate('temp','Dx',Dx,'k',k,'B',B,'N',1,'nullComponent',true,'filters','gabor_4or_32.mat','cc',cc, ...
        'obsVar',sigma,'z_shape',z_shape,'z_scale',z_scale, ...
        'prior','gamma','g_shape',g_shape,'g_scale',g_mean/g_shape,'null_shape',g_shape,'null_scale',g_mean/g_shape);
    models{end+1} = ge2;
    model_names{end+1} = 'gamma2';
    
    ge3 = gestaltCreate('temp','Dx',Dx,'k',k,'B',B,'N',1,'nullComponent',true,'prior','dirichlet','cc',cc, ...
        'filters','gabor_4or_32.mat','obsVar',sigma,'sparsity',0.1,'z_shape',z_shape,'z_scale',z_scale);
    models{end+1} = ge3;
    model_names{end+1} = 'dirichlet';
    
    nModels = size(models,2);    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CREATE STIMULI
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    load('gabor_4or_32.mat')    
    orients = [2 4];
    %orients = [4 2];
    stimuli = {};
    stimuli{1} = 0.1 * randn(Dx,1);
    
%     for l=1:size(locations,1)
%         for o=1:length(orients)
%             %coeffs = zeros(Dx,1);
%             coeffs = randn(Dx,1);
%             for loc=1:size(locations,2)
%                 coeffs((locations(l,loc)-1)*nOrient+orients(o),1) = 5;
%             end
%             stimuli{end+1} = A * coeffs;
%         end
%     end
    
    backgroundZ = 1;
    g_off = 0.01 * ones(k,1);
    g_on = g_off;
    g_on(1,1) = 10;
    X_on = gestaltAncestralSample(ge1,g_on,backgroundZ,false,false);
    X_off = gestaltAncestralSample(ge1,g_off,backgroundZ,false,false);
    stimuli{end+1} = X_on(1,:)';
    stimuli{end+1} = X_off(1,:)';
    
    stimuli{end+1} = 0.1 * randn(Dx,1);
    viewImageSet(stimuli);
    
    exemplar_cells = [(ic_idx-1)*nOrient+2 (ic_idx-1)*nOrient+4 (locations(1)-1)*nOrient+2 (locations(2)-1)*nOrient+4 1];
    exemplar_titles = {'ic1','ic2','stim1','stim2','none'};            
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RUN SAMPLING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    if strcmp(calculation,'transient')
        reset = false;
    elseif strcmp(calculation,'stationary')
        reset = true;
    end
    
    if isempty(appendTo)
        [vsamp,gsamp,zsamp] = gestaltScheduling(stimuli,timings,models,nTrials,0,reset);
    else
        load(appendTo);
    end
        
    if strcmp(calculation,'transient')
        % plot time series
        gdata = permute(gsamp,[1 4 2 3]);
        plotGridSeries(gdata,cumulative_timings,model_names,{},'Model','comp');
        vdata = permute(squeeze(mean(vsamp(:,:,:,:,exemplar_cells),4)),[1 4 2 3]);
        plotGridSeries(vdata,cumulative_timings,model_names,exemplar_titles,'Model','');
        zdata = reshape(zsamp,[nModels 1 nTrials sum(timings)]);
        plotGridSeries(zdata,cumulative_timings,model_names,{'z'},'Model','');
    elseif strcmp(calculation,'stationary')
        burnin = 0;
    end
end

