function [nrf_stim,nrf_allgsamp,cc] = samplingWithLines(nTrials,nSamples,prestimSamp,try_k,k,Dx,filters,v_sampler,randseed,stimFunc,compFunc)
    % reproducing the effect of non-classical stimulation on spike count or 
    % membrane potenital evoked response reliabilities from Haider et al.,
    % Neuron, 2010.
    close all;
    
    if strcmp(randseed,'last')
        load lastrandseed;
    end
    s = RandStream('mt19937ar','Seed',randseed);
    RandStream.setGlobalStream(s);
    randseed = s.Seed;
    save('lastrandseed.mat','randseed');
    
    if strcmp(stimFunc,'reliability')
        createStimulus = @(a,b,c,d) createReliabilityStimuli(a,b,c,d);
    elseif strcmp(stimFunc,'IC')
        createStimulus = @(a,b,c,d,e,f) createICStimuli(a,b,c,d,e,f);
    else
        throw(MException('Gestalt:sampleWithLines:notImplemented','Stimulus generating function %s is not implemented.',stimFunc));
    end
    
    if strcmp(compFunc,'reliability')
        computeFunction = @(a,b,c,d,e) computeReliability(a,b,c,d,e);
    elseif strcmp(stimFunc,'IC')
        computeFunction = @(a,b,c,d,e) computeIC(a,b,c,d,e);
    else
        computeFunction = @(a,b,c,d,e) 1;        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
    % MODEL CREATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % model parameters for generation and sampling
    generating_sigma = 1;
    z_gen = 1;
    g_gen = 10;
        
    sampling_sigma = 1;
    g_scale = 2;
    z_shape = 1;
    z_scale = 0.1;
    sample_z = true;
    fixedZ = 0.1; %in case we don't sample  
    filter_scaling = 10;
    g_init = [zeros(k,1);1];
    
    sampler_verb = 0;
    if strcmp(v_sampler,'mh')
        sampler_verb = 1;
    end
    
    % create model
    if strcmp(filters,'OF')
        filters = sprintf('OF_%d.mat',Dx);
    end
    ge = gestaltCreate('temp','Dx',Dx,'k',k,'B',1,'N',2, ...
        'filters',filters,'obsVar',generating_sigma,'g_scale',g_scale,'z_shape',z_shape,'z_scale',z_scale);
    
    avg_eigenval = mean(abs(eig(ge.A)));
    %ge.A = filter_scaling*(ge.A);
    ge.A = filter_scaling*(1/avg_eigenval)*(ge.A);
    
    % create covariance components that reflect linear shapes
    [~,tmp] = lineImages(100,Dx,k);
    [cc,coeffs] = templates2covariances(tmp,ge.A);
    %viewImageSet(cc);
    cc{k+1} = eye(Dx);
    ge.cc = cc;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
    % STIMULUS CREATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    crf_samples = zeros(try_k,nTrials,nSamples);
    nrf_samples = zeros(try_k,nTrials,nSamples);
    sHandle = figure('Units','normalized','OuterPosition',[0.05 0.4 0.95 0.8]);
    handles = [];
    for c = 1:try_k
        act_coeffs = diag(cc{c});
        ge.obsVar = generating_sigma;
        fprintf('Component %d/%d ',try_k,c);
        [activatedIndices,activatedCoeffs,handles] = getActivatedFilters(cc{c},ge.A,act_coeffs,handles,false);               
        
        [crf_stim,nrf_stim,cell] = createStimulus(ge,activatedIndices,activatedCoeffs,c,z_gen,g_gen);
        
        ge.X(1,:,:) = reshape(crf_stim,1,ge.Dx);
        ge.X(2,:,:) = reshape(nrf_stim,1,ge.Dx);
        
        figure(sHandle);
        clf;        
        nrow = 4;
        ncol = 5;
        
        % plot the stimuli along with cell and gestalt receptive fields        
        subplot(nrow,ncol,3);
        viewImage(ge.A(:,cell),'useMax',true);
        title('Example cell RF');
        subplot(nrow,ncol,4);
        viewImage(tmp{c});
        title('Gestalt RF');
        subplot(nrow,ncol,ncol+3);
        viewImage(crf_stim,'useMax',true);
        title('CRF stimulus');
        xlabel(sprintf('mu = %.3f sigma = %.3f',mean(crf_stim(:)),std(crf_stim(:))));
        subplot(nrow,ncol,ncol+4);        
        viewImage(nrf_stim,'useMax',true);
        title('nCRF stimulus');
        xlabel(sprintf('mu = %.3f sigma = %.3f',mean(nrf_stim(:)),std(nrf_stim(:))));
        %pause;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
        % SAMPLING
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
        crf_gsamp = zeros(nTrials,nSamples);
        nrf_gsamp = zeros(nTrials,nSamples);
        crf_gnullsamp = zeros(nTrials,nSamples);
        nrf_gnullsamp = zeros(nTrials,nSamples);
        crf_zsamp = zeros(nTrials,nSamples);
        nrf_zsamp = zeros(nTrials,nSamples);
        
        crf_activated = zeros(nTrials,length(activatedIndices),nSamples);
        nrf_activated = zeros(nTrials,length(activatedIndices),nSamples);
        crf_allgsamp = zeros(nTrials,ge.k,nSamples);
        nrf_allgsamp = zeros(nTrials,ge.k,nSamples);
        fprintf('%d/%d activated filters ',length(activatedCoeffs),ge.Dv);
        for t = 1:nTrials
            printCounter(t,'stringVal','Trial','maxVal',nTrials,'newLine',true);
            % crf samples
            ge.obsVar = sampling_sigma;
            [cs,~,cz] = gestaltGibbs(ge,1,nSamples,'verbose',0,'vSampler',v_sampler,'contrast',sample_z, ...
                'prestimSamples',prestimSamp,'verbose',sampler_verb,'fixedZ',fixedZ,'initG',g_init);
            crf_samples(c,t,:) = cs(:,ge.k+cell)';
            crf_gsamp(t,:) = cs(:,c);
            crf_gnullsamp(t,:) = cs(:,k+1);
            crf_zsamp(t,:) = cz;
            crf_activated(t,:,:) = cs(:,ge.k+activatedIndices)'; 
            crf_allgsamp(t,:,:) = cs(:,1:ge.k)';
            % nrf samples
            [ns,~,nz] = gestaltGibbs(ge,2,nSamples,'verbose',0,'vSampler',v_sampler,'contrast',sample_z, ...
                'prestimSamples',prestimSamp,'verbose',sampler_verb,'fixedZ',fixedZ,'initG',g_init);
            nrf_samples(c,t,:) = ns(:,ge.k+cell)';
            nrf_gsamp(t,:) = ns(:,c);
            nrf_gnullsamp(t,:) = ns(:,k+1);
            nrf_zsamp(t,:) = nz;
            nrf_activated(t,:,:) = ns(:,ge.k+activatedIndices)';
            nrf_allgsamp(t,:,:) = ns(:,1:ge.k)';
        end                                
        
        figure(sHandle);
        plotPair(squeeze(crf_samples(c,:,:)),squeeze(nrf_samples(c,:,:)),nrow,ncol,1,true,'V',prestimSamp,{});
        plotPair(crf_gsamp,nrf_gsamp,nrow,ncol,1+ncol,true,'G',prestimSamp,{});
        plotPair(crf_gnullsamp,nrf_gnullsamp,nrow,ncol,1+2*ncol,true,'G0',prestimSamp,{});
        plotPair(crf_zsamp,nrf_zsamp,nrow,ncol,1+3*ncol,true,'Z',prestimSamp,{});                        
        
        left_act = squeeze(mean(crf_activated,1));        
        right_act = squeeze(mean(nrf_activated,1));
        labels = {};
        for ii = 1:length(activatedIndices)
            labels{ii} = sprintf('%d',activatedIndices(ii));
        end
        plotPair(left_act,right_act,nrow,ncol,1+2*ncol+2,false,'act-V',prestimSamp,{});                        
        
        left_g = squeeze(mean(crf_allgsamp,1));        
        right_g = squeeze(mean(nrf_allgsamp,1));
        labels = {};
        for ii = 1:ge.k
            labels{ii} = sprintf('%d',ii);
        end
        plotPair(left_g,right_g,nrow,ncol,1+3*ncol+2,false,'all-G',prestimSamp,labels);                        
        
%         if k>1 && c<k
%             pause;
%         end
        
    end
    fprintf('\n');
    computeFunction(crf_samples,nrf_samples,prestimSamp,sHandle,ncol);   
end

function B = upperTriangleValues(A)
    A = A-diag(diag(A));
    B = A(triu(true(size(A))));
    B = B(B~=0);
end

function plotPair(leftData,rightData,plotRows,plotColumns,leftPlotIndex,plotMean,titleString,redLine,labels)
    subplot(plotRows,plotColumns,leftPlotIndex);
    plot(leftData');
    hold on;
    xlim([1,size(leftData,2)]);
    yl1 = ylim();
    if plotMean        
        plot(mean(leftData)','LineWidth',3);
    end
    title(strcat('CRF-',titleString));
    
    subplot(plotRows,plotColumns,leftPlotIndex+1);
    plot(rightData');
    hold on;
    xlim([1,size(rightData,2)]);
    yl2 = ylim();
    if plotMean        
        plot(mean(rightData)','LineWidth',3);
    end
    title(strcat('nCRF-',titleString));
    if ~isempty(labels)
        legend(labels);
    end
    
    subplot(plotRows,plotColumns,leftPlotIndex);
    ylim([min(yl1(1),yl2(1)) max(yl1(2),yl2(2))]);
    plot([redLine;redLine],ylim(),'r-');
    plot([1; size(rightData,2)],[0;0],'k--');
    subplot(plotRows,plotColumns,leftPlotIndex+1);
    ylim([min(yl1(1),yl2(1)) max(yl1(2),yl2(2))]);
    plot([redLine;redLine],ylim(),'r-');
    plot([1; size(rightData,2)],[0;0],'k--');
end

function [crf_stim,nrf_stim,cell] = createICStimuli(ge,activatedIndices,activatedCoeffs,compNum,z_gen,g_gen)
%     [~,maxact] = max(activatedCoeffs);
%     cell = activatedIndices(maxact);
    
    signal_multiplier = 2;

    g = zeros(ge.k,1);
    g(compNum,1) = g_gen;        
    [crf_stim,v] = gestaltAncestralSample(ge,g,z_gen,false,true);
    for i=1:length(activatedIndices)
        crf_stim = crf_stim + signal_multiplier * activatedCoeffs(i) * ge.A(:,activatedIndices(i))';
    end
    
%     
%     v = zeros(ge.Dv,1);
%     v(activatedIndices) = activatedCoeffs;
%     crf_stim = mvnrnd((z*ge.A*v)',ge.obsVar*eye(ge.Dx))'; 
    
     v = v(:);
     [~,maxidx] = max(v(activatedIndices,1));
     cell = activatedIndices(maxidx);
%     
%     v = zeros(ge.Dv,1);
%     v(cell,1) = 1;
%     nrf_stim = mvnrnd((z*ge.A*v)',ge.obsVar*eye(ge.Dx))';     
    
   
    fprintf(' cell idx: %d act: %f ',cell,v(cell,1));
    nrf_stim = crf_stim(:) - (z_gen * v(cell,1) + signal_multiplier*activatedCoeffs(maxidx)) * ge.A(:,cell);
    
end

function [crf_stim,nrf_stim,cell] = createReliabilityStimuli(ge,activatedIndices,activatedCoeffs,compNum)
    % choose a cell     
    [~,maxact] = max(activatedCoeffs);
    cell = activatedIndices(maxact);

    % create a stimulus that lies in the receptive field of the neuron
    z = 1;
    v = zeros(ge.Dv,1);
    v(cell,1) = 1;
    crf_stim = mvnrnd((z*ge.A*v)',ge.obsVar*eye(ge.Dx))';        
%   g = gestaltSamplePriorG(ge,'gamma');
%   [crf_stim,~] = gestaltAncestralSample(ge,g,z,false);        

    % choose a stimulus that also lies in the extraclassical field
    % colinearly with the preferred direction
    g = zeros(ge.k,1);
    g(compNum,1) = 100;        
    [nrf_stim,~] = gestaltAncestralSample(ge,g,z,false);
    nrf_stim = nrf_stim + crf_stim';
%   nrf_stim = createImageStimulus(tmp(c),1);
%   v = zeros(ge.Dv,1);
%   v(cells) = activatedCoeffs;
%   nrf_stim = mvnrnd((z*ge.A*v)',ge.obsVar*eye(ge.Dx))'; 
end

function computeReliability(crf_samples,nrf_samples,prestimSamp,plotHandle,freeSubplots)
    try_k = size(crf_samples,1);
    rel_crf = zeros(try_k,1);
    rel_nrf = zeros(try_k,1);      
    rel_crf_std = zeros(try_k,1);
    rel_nrf_std = zeros(try_k,1);      
    for c = 1:try_k
        corr_crf = corr(squeeze(crf_samples(c,:,:))');
        crf_corrvals = upperTriangleValues(corr_crf);
        rel_crf(c,1) = mean(crf_corrvals);
        rel_crf_std(c,1) = std(crf_corrvals);
        corr_nrf = corr(squeeze(nrf_samples(c,:,:))');
        nrf_corrvals = upperTriangleValues(corr_nrf);
        rel_nrf(c,1) = mean(nrf_corrvals);
        rel_nrf_std(c,1) = std(nrf_corrvals);
    end
    fprintf('\n');
    mean_crf = mean(rel_crf);
    std_crf = mean(rel_crf_std);
    mean_nrf = mean(rel_nrf);
    std_nrf = mean(rel_nrf_std);
    
    % plot results
    figure(plotHandle);
    subplot(3,4,freeSubplots(1));
    barwitherr([std_crf std_nrf],[mean_crf mean_nrf]);
    ymin = min([0,mean_crf,mean_nrf]);
    if ymin < 0
        ymin = min([mean_crf - std_crf - 0.1,mean_nrf - std_nrf - 0.1]);
    end
    ymax = max([0.35,mean_crf + std_crf + 0.1,mean_nrf + std_nrf + 0.1]);
    ylim([ymin ymax]);
    set(gca,'XTickLabel',{'CRF','nCRF'});
    title(sprintf('K = %d',try_k));
    ylabel('reliability');
end

function computeIC(crf_samples,nrf_samples,prestimSamp,plotHandle,freeSubplots)
    real = squeeze(mean(crf_samples,1));
    ic = squeeze(mean(nrf_samples,1));
    real_mean = mean(real,1);
    real_std = std(real,1);
    ic_mean = mean(ic,1);
    ic_std = std(ic,1);
    
    first = prestimSamp + 1;
    last = prestimSamp + 3;
    
    figure(plotHandle);
    subplot(4,5,freeSubplots(1));
    barwitherr([real_std(first:last)' ic_std(first:last)'],[real_mean(first:last)' ic_mean(first:last)']);
    
end