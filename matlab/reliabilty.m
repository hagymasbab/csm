function reliabilty(nTrials,nSamples,k,Dx,filters)
    % reproducing the effect of non-classical stimulation on spike count or 
    % membrane potenital evoked response reliabilities from Haider et al.,
    % Neuron, 2010.
    close all;
    
    % model parameters for generation and sampling
    generating_sigma = 0.001;
    sampling_sigma = 1;
    g_scale = 2;
    
    % create model
    if strcmp(filters,'OF')
        filters = sprintf('OF_%d.mat',Dx);
    end
    ge = gestaltCreate('temp','Dx',Dx,'k',k,'B',1,'filters',filters,'obsVar',generating_sigma,'N',2,'g_scale',g_scale);
    
    % create covariance components that reflect linear shapes
    [~,tmp] = lineImages(100,Dx,k);
    [cc,coeffs] = templates2covariances(tmp,ge.A);
    %viewImageSet(cc);
    cc{k+1} = eye(Dx);
    ge.cc = cc;
            
    rel_crf = zeros(k,1);
    rel_nrf = zeros(k,1);      
    sHandle = figure('Units','normalized','OuterPosition',[0.1 0.4 0.8 0.4]);
    pfHandle = figure('Units','normalized','OuterPosition',[0.1 0.1 0.4 0.8]);
    nfHandle = figure('Units','normalized','OuterPosition',[0.6 0.1 0.4 0.8]);
    for c = 1:k
        ge.obsVar = generating_sigma;
        fprintf('Component %d/%d ',k,c);
        
        % select the neurons that are activated by the component
        actFilters = covariance2template(cc{k},ge.A);
        cell = find(actFilters);
        activatedFilters = ge.A(:,cell)';
        activatedCoeffs = coeffs(c,cell)';
        positiveFilters = activatedFilters(activatedCoeffs>0,:);
        negativeFilters = activatedFilters(activatedCoeffs<0,:);
        positiveCoeffs = activatedCoeffs(activatedCoeffs>0,1);
        negativeCoeffs = activatedCoeffs(activatedCoeffs<0,1);
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
        numfilter = length(cell);
        
        % choose a cell     
        [~,maxact] = max(activatedCoeffs);
        cell = cell(maxact);
        
        % create a stimulus that lies in the receptive field of the neuron
        v = zeros(ge.Dv,1);
        v(cell,1) = 1;
        crf_stim = mvnrnd((ge.A*v)',ge.obsVar*eye(ge.Dx))';        
        ge.X(1,:,:) = reshape(crf_stim,1,ge.Dx);
        
        % choose a stimulus that also lies in the extraclassical field
        % colinearly with the preferred direction
        g = zeros(ge.k,1);
        g(c,1) = 1;        
        [nrf_stim,~] = gestaltAncestralSample(ge,g,1,false);
%         nrf_stim = createImageStimulus(tmp(c),1);
        ge.X(2,:,:) = reshape(nrf_stim,1,ge.Dx);
        
        crf_samples = zeros(nTrials,nSamples);
        nrf_samples = zeros(nTrials,nSamples);
        crf_gsamp = zeros(nTrials,nSamples);
        nrf_gsamp = zeros(nTrials,nSamples);
        fprintf('%d/%d activated filters ',numfilter,ge.Dv);
        for t = 1:nTrials
            %fprintf('\nTrial %d/%d\n',nTrials,t);
            printCounter(t,'stringVal','Trial','maxVal',nTrials,'newLine',true);
            % crf samples
            ge.obsVar = sampling_sigma;
            cs = gestaltGibbs(ge,1,nSamples,'verbose',0);
            crf_samples(t,:) = cs(:,ge.k+cell)';
            crf_gsamp(t,:) = cs(:,c);
            % nrf samples
            ns = gestaltGibbs(ge,2,nSamples,'verbose',0);
            nrf_samples(t,:) = ns(:,ge.k+cell)';
            nrf_gsamp(t,:) = ns(:,c);
        end        
        
        figure(sHandle);
        clf;
        subplot(2,4,1);
        plot(crf_samples');
        xlim([1,nSamples]);
        subplot(2,4,2);
        plot(nrf_samples');
        xlim([1,nSamples]);
        subplot(2,4,5);
        plot(crf_gsamp');
        xlim([1,nSamples]);
        hold on;
        plot(mean(crf_gsamp)','LineWidth',3);
        subplot(2,4,6);
        plot(nrf_gsamp');
        xlim([1,nSamples]);
        hold on;
        plot(mean(nrf_gsamp)','LineWidth',3);
        
        % plot the stimuli along with cell and gestalt receptive fields        
        subplot(2,4,3);
        viewImage(ge.A(:,cell),'useMax',true);
        title('Example cell RF');
        subplot(2,4,4);
        viewImage(tmp{c});
        title('Gestalt RF');
        subplot(2,4,7);
        viewImage(crf_stim,'useMax',true);
        title('CRF stimulus');
        subplot(2,4,8);
        viewImage(nrf_stim,'useMax',true);
        title('nCRF stimulus');
        
        if k>1 && c<k
            pause;
        end
        
        corr_crf = corr(crf_samples');
        %size(corr_crf)
        rel_crf(c,1) = mean(upperTriangleValues(corr_crf));
        corr_nrf = corr(nrf_samples');
        rel_nrf(c,1) = mean(upperTriangleValues(corr_nrf));
    end
    fprintf('\n');
    % plot results
    figure();
    barwitherr([std(rel_crf) std(rel_nrf)],[mean(rel_crf) mean(rel_nrf)]);
    set(gca,'XTickLabel',{'CRF','nCRF'});
    title(sprintf('N = %d',k));
    ylabel('reliability');
   
end

function B = upperTriangleValues(A)
    A = A-diag(diag(A));
    B = A(triu(true(size(A))));
    B = B(B~=0);
end