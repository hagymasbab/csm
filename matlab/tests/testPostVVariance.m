function testPostVVariance(Dv,k,sampleSizes,nTrials,loadVars,plotHist)
    % for all sample sizes
    % calculate nTrials posterior covariance matrices with different
    % samples
    % calculate average variance of matrix elemnts (upper triangle)
    
    ge = gestaltCreate('temp','Dx',Dv,'k',k,'filters','gabor_4or','obsVar',0.5,'N',1, ...
        'g_shape',2,'g_scale',2,'z_shape',2,'z_scale',2,'N',1,'generateComponents',true,'generateData',true);
    
    x = reshape(ge.X(1,1,:),ge.Dx,1);
    
    if loadVars
        load('save_testcorrvar.mat');
    else
        variabilities = cell(1,length(sampleSizes));
        for i=1:length(sampleSizes)
            corrmats = zeros(nTrials,ge.Dv,ge.Dv);
            for t = 1:nTrials
                covc = posteriorCovariances(x,ge,sampleSizes(i),'shuffle',false);
                corrmats(t,:,:) = corrcov(covc);
            end
            stds = squeeze(std(corrmats,0,1));
            variabilities{i} = upperTriangleValues(stds);
        end
        save('bin/save_testcorrvar.mat','variabilities');
    end
    
    if plotHist        
        vm = cell2mat(variabilities);
        maxvar = max(vm(:));
        histbins = linspace(0,maxvar,22);
        histbins = histbins(2:end-1);
        maxcount = 0;
        for i=1:length(sampleSizes)
            subplot(1,length(sampleSizes),i);
            counts = hist(variabilities{i},histbins);            
            if max(counts) > maxcount 
                maxcount = max(counts);
            end
            hist(variabilities{i},histbins);
            xlim([0 maxvar])
            xlabel(sprintf('# of samples = %d',sampleSizes(i)),'FontSize',16)
        end
        for i=1:length(sampleSizes)
            subplot(1,length(sampleSizes),i);
            ylim([0 maxcount*1.1]);
        end
    end
end