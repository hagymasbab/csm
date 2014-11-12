function reliabilty(nTrials,nSamples,k,Dx)
    % reproducing the effect of non-classical stimulation on spike count or 
    % membrane potenital evoked response reliabilities from Haider et al.,
    % Neuron, 2010.
    close all;
    
    % create model
    ge = gestaltCreate('temp','Dx',Dx,'k',k,'B',1,'filters','eye','obsVar',0.001,'N',2);
    
    % create covariance components that reflect linear shapes
    [~,tmp] = lineImages(100,Dx,k);
    [cc,~] = templates2covariances(tmp,ge.A);
    %viewImageSet(cc);
    cc{k+1} = eye(Dx);
    ge.cc = cc;
            
    rel_crf = zeros(k,1);
    rel_nrf = zeros(k,1);      
    figure('Units','normalized','OuterPosition',[0.1 0.6 0.8 0.4]);
    for c = 1:k
        ge.obsVar = 0.01;
        fprintf('Component %d/%d',k,c);
        % select the neurons that are activated by the component
        actFilters = gestaltExtractTemplate(ge,k);
        % choose a cell 
        cell = find(actFilters);
        cell = cell(1);
        % create a stimulus that lies in the receptive field of the neuron
        v = zeros(ge.Dv,1);
        v(cell,1) = 1;
        crf_stim = mvnrnd((ge.A*v)',ge.obsVar*eye(ge.Dx))';        
        ge.X(1,:,:) = reshape(crf_stim,1,ge.Dx);
        % choose a stimulus that also lies in the extraclassical field
        % colinearly with the preferred direction
%         g = zeros(ge.k,1);
%         g(c,1) = 1;        
%         [nrf_stim,~] = gestaltAncestralSample(ge,g,1,false);
        nrf_stim = createImageStimulus(tmp(c),1);
        ge.X(2,:,:) = reshape(nrf_stim,1,ge.Dx);
        
        crf_samples = zeros(nTrials,nSamples);
        nrf_samples = zeros(nTrials,nSamples);
        for t = 1:nTrials
            fprintf('\nTrial %d/%d\n',nTrials,t);
            % crf samples
            ge.obsVar = 1;
            cs = gestaltGibbs(ge,1,nSamples,'verbose',1);
            crf_samples(t,:) = cs(:,ge.k+cell)';
            % nrf samples
            ns = gestaltGibbs(ge,2,nSamples,'verbose',1);
            nrf_samples(t,:) = ns(:,ge.k+cell)';
        end
                        
        subplot(1,2,1);
        plot(crf_samples');
        xlim([1,nSamples]);
        subplot(1,2,2);
        plot(nrf_samples');
        xlim([1,nSamples]);
        
        corr_crf = corr(crf_samples');
        %size(corr_crf)
        rel_crf(c,1) = mean(upperTriangleValues(corr_crf));
        corr_nrf = corr(nrf_samples');
        rel_nrf(c,1) = mean(upperTriangleValues(corr_nrf));
    end
    fprintf('\n');
    % plot results
    figure();
    subplot(2,3,1);
    barwitherr([std(rel_crf) std(rel_nrf)],[mean(rel_crf) mean(rel_nrf)]);
    set(gca,'XTickLabel',{'CRF','CRF+nCRF'});
    title(sprintf('N = %d',k));
    ylabel('reliability');
    % plot example stimuli along with cell and gestalt receptive fields
    % we will plot the last one computed
    subplot(2,3,2);
    viewImage(ge.A(:,cell));
    title('Example cell RF');
    subplot(2,3,3);
    viewImage(tmp{k});
    title('Gestalt RF');
    subplot(2,3,4);
    viewImage(crf_stim);
    title('CRF stimulus');
    subplot(2,3,5);
    viewImage(nrf_stim);
    title('nCRF stimulus');
end

function B = upperTriangleValues(A)
    A = A-diag(diag(A));
    B = A(triu(true(size(A))));
    B = B(B~=0);
end