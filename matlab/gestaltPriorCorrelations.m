function gestaltPriorCorrelations(nTrials)
    Dx = 1024;    
    z_shape = 1;
    z_scale = 0.1;
    sigma = 1;
    g_mean = 1;

    % TODO create stimuli
    stimuli = {};
    stimuli{1} = 1;
    
    
    % define timings
    timings = [3,5,5,3];
    cumulative_timings = cumsum(timings);
    
    % TODO create covariance components
    filterlists = [];
    cc = filterList2Components(filterlists,true,Dx);
    k = size(cc,2); 

    models = {};
    model_names = {};
    
    g_shape = 1;
    ge1 = gestaltCreate('temp','Dx',Dx,'k',k,'B',1,'N',1,'nullComponent',true,'filters','gabor_4or_32.mat','cc',cc, ...
        'obsVar',sigma,'z_shape',z_shape,'z_scale',z_scale, ...
        'prior','gamma','g_shape',g_shape,'g_scale',g_mean/g_shape,'null_shape',g_shape,'null_scale',g_mean/g_shape);
    models{end+1} = ge1;
    model_names{end+1} = 'gamma_1';
    
    g_shape = 2;
    ge2 = gestaltCreate('temp','Dx',Dx,'k',k,'B',1,'N',1,'nullComponent',true,'filters','gabor_4or_32.mat','cc',cc, ...
        'obsVar',sigma,'z_shape',z_shape,'z_scale',z_scale, ...
        'prior','gamma','g_shape',g_shape,'g_scale',g_mean/g_shape,'null_shape',g_shape,'null_scale',g_mean/g_shape);
    models{end+1} = ge2;
    model_names{end+1} = 'gamma_2';
    
    ge3 = gestaltCreate('temp','Dx',Dx,'k',k,'B',1,'N',1,'nullComponent',true,'prior','dirichlet','cc',cc, ...
        'filters','gabor_4or_32.mat','obsVar',sigma,'sparsity',0.1,'z_shape',z_shape,'z_scale',z_scale);
    models{end+1} = ge3;
    model_names{end+1} = 'dirichlet';
    
    nModels = size(models,2);
    
    [vsamp,gsamp,zsamp] = gestaltScheduling(stimuli,timings,models,nTrials);
    
    % plot g time series
    for m = 1:nModels
        for kk = 1:k
            actdata = gsamp(m,:,:,kk);
            subplot(k,nModels,(kk-1)*nModels+m);
            plot(actdata');
            hold on;
            xlim([1,size(actdata,2)]);
            ylim([0 1]);
            plot(mean(actdata)','LineWidth',3);            
            title(sprintf('Model %s comp %d',model_names{m},kk),'FontSize',16);
            for t=1:length(timings)
                act_t = cumulative_timings(t);
                plot([act_t;act_t],ylim(),'r-');
            end
        end
    end
    
end