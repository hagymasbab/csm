function gestaltPriorCorrelations(nTrials,timings,appendTo)
    close all;
    Dx = 1024;    
    
    nOrient = 4; % this is true for gabor_4or_32.mat only
    shift = nOrient/2; 
    rfsinarow = sqrt(Dx)/shift;    
    
    z_shape = 1;
    z_scale = 0.1;
    sigma = 0.1;
    g_mean = 0.1;

    % create stimuli
    % backgroundZ = 1;
    load('gabor_4or_32.mat')
    ic_idx = floor(Dx / (nOrient*2) + rfsinarow/2);
    sh1 = 3; % single shift 
    sh2 = 5; % double shift
    locations = [ic_idx-(sh2*rfsinarow)+sh2 ic_idx-(sh1*rfsinarow)+sh1 ic_idx+(sh1*rfsinarow)-sh1 ic_idx+(sh2*rfsinarow)-sh2];
    orients = [2 4];
    %orients = [4 2];
    stimuli = {};
    stimuli{1} = 0.1 * randn(Dx,1);
    for l=1:size(locations,1)
        for o=1:length(orients)
            %coeffs = zeros(Dx,1);
            coeffs = randn(Dx,1);
            for loc=1:size(locations,2)
                coeffs((locations(l,loc)-1)*nOrient+orients(o),1) = 5;
            end
            stimuli{end+1} = A * coeffs + 0.1 * randn(Dx,1);
        end
    end
    stimuli{end+1} = 0.1 * randn(Dx,1);
    
    exemplar_cells = [(ic_idx-1)*nOrient+2 (ic_idx-1)*nOrient+4 (locations(1)-1)*nOrient+2 (locations(2)-1)*nOrient+4 1];
    exemplar_titles = {'ic1','ic2','stim1','stim2','none'};
    
    %viewImageSet(stimuli,'max',false);
    
    % define timings
    % timings = [20,80,80,20];
    cumulative_timings = cumsum(timings);
    
    % create covariance components
    filterlists = [ ([ic_idx locations] - 1) .* nOrient + 2; ([ic_idx locations] - 10) .* nOrient + 3 ];
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
    models{end+1} = ge3;
    model_names{end+1} = 'dirichlet';
    
    nModels = size(models,2);        
    
    if isempty(appendTo)
        [vsamp,gsamp,zsamp] = gestaltScheduling(stimuli,timings,models,nTrials);
    else
        load(appendTo);
    end

    
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
                plot([act_t;act_t],ylim(),'r-','LineWidth',3);
            end
        end
    end
    
    % plot v time series
    figure();
    for m = 1:nModels
        for c = 1:length(exemplar_cells)
            actdata = squeeze(vsamp(m,:,:,exemplar_cells(c)));
            subplot(length(exemplar_cells),nModels,(c-1)*nModels+m);
            plot(actdata');
            hold on;
            xlim([1 sum(timings)]);
            ylim([-4 8]);
            plot(mean(actdata)','LineWidth',3);            
            title(sprintf('Model %s %s',model_names{m},exemplar_titles{c}),'FontSize',16);
            for t=1:length(timings)-1
                act_t = cumulative_timings(t);
                plot([act_t;act_t],ylim(),'r-','LineWidth',3);
            end
        end
    end
    
    % plot z time series
    figure();
    for m = 1:nModels        
        actdata = squeeze(zsamp(m,:,:));
        subplot(1,nModels,m);
        plot(actdata');
        hold on;
        xlim([1 sum(timings)]);
        ylim([0 1]);
        plot(mean(actdata)','LineWidth',3);            
        title(sprintf('Model %s',model_names{m}),'FontSize',16);
        for t=1:length(timings)-1
            act_t = cumulative_timings(t);
            plot([act_t;act_t],ylim(),'r-','LineWidth',3);
        end        
    end
    
end