function [mean_crf,std_crf,mean_nrf,std_nrf] = reliabilty(nTrials,nSamples,try_k,k,Dx,filters)
    % reproducing the effect of non-classical stimulation on spike count or 
    % membrane potenital evoked response reliabilities from Haider et al.,
    % Neuron, 2010.
    close all;
    
    % model parameters for generation and sampling
    generating_sigma = 0.001;
    sampling_sigma = 1;
    g_scale = 2;
    z_shape = 1;
    z_scale = 0.1;
    v_sampler = 'direct';
    sample_z = true;
    
    % create model
    if strcmp(filters,'OF')
        filters = sprintf('OF_%d.mat',Dx);
    end
    ge = gestaltCreate('temp','Dx',Dx,'k',k,'B',1,'N',2, ...
        'filters',filters,'obsVar',generating_sigma,'g_scale',g_scale,'z_shape',z_shape,'z_scale',z_scale);
    
    % create covariance components that reflect linear shapes
    [~,tmp] = lineImages(100,Dx,k);
    [cc,coeffs] = templates2covariances(tmp,ge.A);
    %viewImageSet(cc);
    cc{k+1} = eye(Dx);
    ge.cc = cc;
            
    rel_crf = zeros(try_k,1);
    rel_nrf = zeros(try_k,1);      
    rel_crf_std = zeros(try_k,1);
    rel_nrf_std = zeros(try_k,1);      
    sHandle = figure('Units','normalized','OuterPosition',[0.1 0.4 0.8 0.8]);
    pfHandle = figure('Units','normalized','OuterPosition',[0.05 0.1 0.4 0.8]);
    nfHandle = figure('Units','normalized','OuterPosition',[0.55 0.1 0.4 0.8]);
    for c = 1:try_k
        ge.obsVar = generating_sigma;
        fprintf('Component %d/%d ',try_k,c);
        
        % select the neurons that are activated by the component
        actFilters = covariance2template(cc{k},ge.A);
        cells = find(actFilters);
        activatedFilters = ge.A(:,cells)';
        activatedCoeffs = coeffs(c,cells)';
        positiveFilters = activatedFilters(activatedCoeffs>0,:);
        negativeFilters = activatedFilters(activatedCoeffs<0,:);
        [positiveCoeffs,posperm] = sort(activatedCoeffs(activatedCoeffs>0,1),1,'descend');
        [negativeCoeffs,negperm] = sort(activatedCoeffs(activatedCoeffs<0,1),1,'ascend');        
        positiveFilters = positiveFilters(posperm,:);
        negativeFilters = negativeFilters(negperm,:);
        positiveTitles = {};
        negativeTitles = {};
        for i=1:length(positiveCoeffs)
            positiveTitles{i} = num2str(positiveCoeffs(i,1));
        end
        for i=1:length(negativeCoeffs)
            negativeTitles{i} = num2str(negativeCoeffs(i,1));
        end        
        figure(pfHandle);        
        viewImageSet(positiveFilters,'titles',positiveTitles);
        figure(nfHandle);
        viewImageSet(negativeFilters,'titles',negativeTitles);
        numfilter = length(cells);
        
        % choose a cell     
        [~,maxact] = max(activatedCoeffs);
        cell = cells(maxact);
        
        % create a stimulus that lies in the receptive field of the neuron
        z = 1;
        v = zeros(ge.Dv,1);
        v(cell,1) = 1;
        crf_stim = mvnrnd((z*ge.A*v)',ge.obsVar*eye(ge.Dx))';        
%         g = gestaltSamplePriorG(ge,'gamma');
%         [crf_stim,~] = gestaltAncestralSample(ge,g,z,false);        
        ge.X(1,:,:) = reshape(crf_stim,1,ge.Dx);
        
        % choose a stimulus that also lies in the extraclassical field
        % colinearly with the preferred direction
        g = zeros(ge.k,1);
        g(c,1) = 100;        
        [nrf_stim,~] = gestaltAncestralSample(ge,g,z,false);
        nrf_stim = nrf_stim + crf_stim';
%         nrf_stim = createImageStimulus(tmp(c),1);
%         v = zeros(ge.Dv,1);
%         v(cells) = activatedCoeffs;
%         nrf_stim = mvnrnd((z*ge.A*v)',ge.obsVar*eye(ge.Dx))'; 
        ge.X(2,:,:) = reshape(nrf_stim,1,ge.Dx);
        
        figure(sHandle);
        clf;        
        nrow = 4;
        ncol = 4;
        
        % plot the stimuli along with cell and gestalt receptive fields        
        subplot(nrow,ncol,3);
        viewImage(ge.A(:,cell),'useMax',true);
        title('Example cell RF');
        subplot(nrow,ncol,4);
        viewImage(tmp{c});
        title('Gestalt RF');
        subplot(nrow,ncol,7);
        viewImage(crf_stim,'useMax',true);
        title('CRF stimulus');
        xlabel(sprintf('mu = %.3f sigma = %.3f',mean(crf_stim(:)),std(crf_stim(:))));
        subplot(nrow,ncol,8);        
        viewImage(nrf_stim,'useMax',true);
        title('nCRF stimulus');
        xlabel(sprintf('mu = %.3f sigma = %.3f',mean(nrf_stim(:)),std(nrf_stim(:))));
        
        crf_samples = zeros(nTrials,nSamples);
        nrf_samples = zeros(nTrials,nSamples);
        crf_gsamp = zeros(nTrials,nSamples);
        nrf_gsamp = zeros(nTrials,nSamples);
        crf_gnullsamp = zeros(nTrials,nSamples);
        nrf_gnullsamp = zeros(nTrials,nSamples);
        crf_zsamp = zeros(nTrials,nSamples);
        nrf_zsamp = zeros(nTrials,nSamples);
        fprintf('%d/%d activated filters ',numfilter,ge.Dv);
        for t = 1:nTrials
            %fprintf('\nTrial %d/%d\n',nTrials,t);
            printCounter(t,'stringVal','Trial','maxVal',nTrials,'newLine',true);
            % crf samples
            ge.obsVar = sampling_sigma;
            [cs,~,cz] = gestaltGibbs(ge,1,nSamples,'verbose',0,'vSampler',v_sampler,'contrast',sample_z);
            crf_samples(t,:) = cs(:,ge.k+cell)';
            crf_gsamp(t,:) = cs(:,c);
            crf_gnullsamp(t,:) = cs(:,k+1);
            crf_zsamp(t,:) = cz;
            % nrf samples
            [ns,~,nz] = gestaltGibbs(ge,2,nSamples,'verbose',0,'vSampler',v_sampler,'contrast',sample_z);
            nrf_samples(t,:) = ns(:,ge.k+cell)';
            nrf_gsamp(t,:) = ns(:,c);
            nrf_gnullsamp(t,:) = ns(:,k+1);
            nrf_zsamp(t,:) = nz;
        end                        
        
        plotPair(crf_samples,nrf_samples,nrow,ncol,1,false,'V');
        plotPair(crf_gsamp,nrf_gsamp,nrow,ncol,1+ncol,true,'G');
        plotPair(crf_gnullsamp,nrf_gnullsamp,nrow,ncol,1+2*ncol,true,'G0');
        plotPair(crf_zsamp,nrf_zsamp,nrow,ncol,1+3*ncol,true,'Z');                        
        
%         if k>1 && c<k
%             pause;
%         end
        
        corr_crf = corr(crf_samples');
        crf_corrvals = upperTriangleValues(corr_crf);
        rel_crf(c,1) = mean(crf_corrvals);
        rel_crf_std(c,1) = std(crf_corrvals);
        corr_nrf = corr(nrf_samples');
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
    figure(sHandle);
    subplot(3,4,11);
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

function B = upperTriangleValues(A)
    A = A-diag(diag(A));
    B = A(triu(true(size(A))));
    B = B(B~=0);
end

function plotPair(leftData,rightData,plotRows,plotColumns,leftPlotIndex,plotMean,titleString)
    subplot(plotRows,plotColumns,leftPlotIndex);
    plot(leftData');
    xlim([1,size(leftData,2)]);
    yl1 = ylim();
    if plotMean
        hold on;
        plot(mean(leftData)','LineWidth',3);
    end
    title(strcat('CRF-',titleString));
    
    subplot(plotRows,plotColumns,leftPlotIndex+1);
    plot(rightData');
    xlim([1,size(rightData,2)]);
    yl2 = ylim();
    if plotMean
        hold on;
        plot(mean(rightData)','LineWidth',3);
    end
    title(strcat('nCRF-',titleString));
    
    subplot(plotRows,plotColumns,leftPlotIndex);
    ylim([min(yl1(1),yl2(1)) max(yl1(2),yl2(2))]);
    subplot(plotRows,plotColumns,leftPlotIndex+1);
    ylim([min(yl1(1),yl2(1)) max(yl1(2),yl2(2))]);
end