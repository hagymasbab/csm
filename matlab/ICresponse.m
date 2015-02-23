function ICresponse(nTrials,nSamples,loadSamples,randseed)
    close all;        
    setrandseed(randseed);
    
    Dx = 1024;
    k = 2;
    nOrient = 4; 
    shift = nOrient/2; 
    rfsinarow = sqrt(Dx)/shift;  
    ic_idx = floor(Dx / (nOrient*2) + rfsinarow/2);
    sh1 = 2; % single shift 
    sh2 = 4;        
        
    central_cell = 8*64+25 + 1;
    upper_cell = 12*64+9 + 1;
    omitted_cell = 10*64+17 + 1;
    lower_cell = 4*64+41 + 1;
    other_cell = 6*64+33 + 1;
    gestalts = [central_cell upper_cell omitted_cell;central_cell lower_cell other_cell];
    cc = filterList2Components(gestalts,false,Dx);
    %viewImageSet(cc);    
    
    g_shape = 1;
    g_scale = 1;
    z_shape = 1;
    z_scale = 0.1;
    sigma = 1;    
    
    ge = gestaltCreate('temp','Dx',Dx,'k',k,'B',1,'N',1,'nullComponent',false,'filters','gabor_4or','cc',cc, ...
            'obsVar',sigma,'z_shape',z_shape,'z_scale',z_scale,'g_shape',g_shape,'g_scale',g_scale);
    
    stimuli = cell(1,2);
    load(sprintf('filters_gabor_4or_%d.mat',Dx));
    mixer = 10;
    coeffs = zeros(Dx,1);    
    coeffs(central_cell,1) = mixer;
    coeffs(upper_cell,1) = mixer;
    %coeffs(omitted_cell,1) = mixer;
    stimuli{1} = A * coeffs;
    coeffs = zeros(Dx,1);    
    coeffs(central_cell,1) = mixer;
    coeffs(lower_cell,1) = mixer;
    %coeffs(other_cell,1) = mixer;
    stimuli{2} = A * coeffs;
    viewImageSet(stimuli);
    
    timings = [nSamples nSamples];
        
    if loadSamples
        load('icsamples.mat');        
    else
        [vsamp,gsamp,zsamp] = gestaltScheduling(stimuli,timings,{ge},nTrials,0,false,'gibbs',false);
        save('icsamples.mat','vsamp','gsamp','zsamp');
    end
    
    exemplar_cells = [upper_cell omitted_cell lower_cell other_cell];
    exemplar_titles = {'stim1','ic1','stim2','ic2'};
    
    cumulative_timings = cumsum(timings);
    gdata = permute(gsamp(:,:,:,1:2),[1 4 2 3]);
    plotGridSeries(gdata,cumulative_timings,{},{'1','2'},'','comp');
    vdata = permute(reshape(mean(vsamp(:,:,:,:,exemplar_cells),4),1,nTrials,cumulative_timings(end),length(exemplar_cells)),[1 4 2 3]);
    plotGridSeries(vdata,cumulative_timings,{''},exemplar_titles,'','');
end