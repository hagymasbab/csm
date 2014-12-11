function gestaltPriorCorrelations(nTrials)
    close all;
    Dx = 1024;    
    
    nOrient = 4; % this is true for gabor_4or_32.mat only
    shift = nOrient/2; 
    rfsinarow = sqrt(Dx)/shift;    
    
    z_shape = 1;
    z_scale = 0.1;
    sigma = 1;
    g_mean = 1;

    % create stimuli
    % backgroundZ = 1;
    load('gabor_4or_32.mat')
    ic_idx = floor(Dx / (nOrient*2) + rfsinarow/2);
    sh1 = 3; % single shift 
    sh2 = 5; % double shift
    locations = [ic_idx-(sh2*rfsinarow)+sh2 ic_idx-(sh1*rfsinarow)+sh1 ic_idx+(sh1*rfsinarow)-sh1 ic_idx+(sh2*rfsinarow)-sh2];
    orients = [2 4];
    stimuli = {};
    stimuli{1} = 0.1 * randn(Dx,1);
    for l=1:size(locations,1)
        for o=1:length(orients)
            coeffs = zeros(Dx,1);
            for loc=1:size(locations,2)
                coeffs((locations(l,loc)-1)*nOrient+orients(o),1) = 1;
            end
            stimuli{end+1} = A * coeffs;
        end
    end
    stimuli{end+1} = 0.1 * randn(Dx,1);
    viewImageSet(stimuli,'max',false);
    
    % define timings
    timings = [3,10,10,3];
    cumulative_timings = cumsum(timings);
    
    % create covariance components
    filterlists = ([ic_idx locations] - 1) .* nOrient + 2;
    cc = filterList2Components(filterlists,true,Dx);
    k = size(cc,2); 
    %viewImageSet(cc);

    models = {};
    model_names = {};
    
    g_shape = 1;
    ge1 = gestaltCreate('temp','Dx',Dx,'k',k,'B',1,'N',1,'nullComponent',true,'filters','gabor_4or_32.mat','cc',cc, ...
        'obsVar',sigma,'z_shape',z_shape,'z_scale',z_scale, ...
        'prior','gamma','g_shape',g_shape,'g_scale',g_mean/g_shape,'null_shape',g_shape,'null_scale',g_mean/g_shape);
    models{end+1} = ge1;
    model_names{end+1} = 'gamma1';
    
    g_shape = 2;
    ge2 = gestaltCreate('temp','Dx',Dx,'k',k,'B',1,'N',1,'nullComponent',true,'filters','gabor_4or_32.mat','cc',cc, ...
        'obsVar',sigma,'z_shape',z_shape,'z_scale',z_scale, ...
        'prior','gamma','g_shape',g_shape,'g_scale',g_mean/g_shape,'null_shape',g_shape,'null_scale',g_mean/g_shape);
    models{end+1} = ge2;
    model_names{end+1} = 'gamma2';
    
    ge3 = gestaltCreate('temp','Dx',Dx,'k',k,'B',1,'N',1,'nullComponent',true,'prior','dirichlet','cc',cc, ...
        'filters','gabor_4or_32.mat','obsVar',sigma,'sparsity',0.1,'z_shape',z_shape,'z_scale',z_scale);
    %models{end+1} = ge3;
    model_names{end+1} = 'dirichlet';
    
    nModels = size(models,2);
    
    [vsamp,gsamp,zsamp] = gestaltScheduling(stimuli,timings,models,nTrials);
    
    % plot g time series
    figure();
    for m = 1:nModels
        for kk = 1:k
            actdata = squeeze(gsamp(m,:,:,kk));
            subplot(k,nModels,(kk-1)*nModels+m);
            plot(actdata');
            hold on;
            xlim([1 sum(timings)]);
            ylim([0 1]);
            plot(mean(actdata)','LineWidth',3);            
            title(sprintf('Model %s comp %d',model_names{m},kk),'FontSize',16);
            for t=1:length(timings)-1
                act_t = cumulative_timings(t);
                plot([act_t;act_t],ylim(),'r-');
            end
        end
    end
    
end