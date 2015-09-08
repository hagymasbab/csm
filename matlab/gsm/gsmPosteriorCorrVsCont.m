function gsmPosteriorCorrVsCont(A,C,sigma_x,zs,x_base)
    setrandseed(1);
    nZ = length(zs);
    iC = stableInverse(C);
    sATA = (1/sigma_x^2) * (A' * A);
    
%     if isempty(x_base)        
        
%         ymax = 0;
%         for i=1:nZ
%             actCov = stableInverse(iC + zs(i)^2 * sATA);
%             actCorr = upperTriangleValues(corrcov(actCov));
%             subplot(1,nZ,i);
%             hist(actCorr,linspace(-1,1,100));
%             yl = ylim();
%             if yl(2) > ymax
%                 ymax = yl(2);
%             end
%         end
%         for i=1:nZ
%             subplot(1,nZ,i);
%             xlim([-1 1]);
%             ylim([0 ymax*1.05]);
%         end
%     else

    if ischar(x_base) 
        if strcmp(x_base,'random')
            x_base = randn(size(A,1),1);
        elseif strcmp(x_base,'gen')
            x_base = gsmGenerate(1,C,A,sigma_x,2,2)';
        else
            error('e');
        end
    end
    pmax = 0;
    rmax = 0;
    smax = 0;
    nSamp = 20000;
    for i=1:nZ
        printCounter(i,'maxVal',nZ,'stringVal','Contrast');

        if isempty(x_base)        
            C_post = stableInverse(iC + zs(i)^2 * sATA);
            samples = mvnrnd(zeros(nSamp,size(A,2)),C_post);

        else
            % construct x        
            x_act = zs(i) * x_base;
            % compute posterior
            [mu_post,C_post,~,~,z_post_dens,component_mus,component_Cs] = gsmPosteriorV(x_act,A,C,sigma_x,2,2,10);

            if any(isnan(z_post_dens))
                fprintf('NaN is Z dens\n');
                continue
            end
            % sample back from mixture
            samples = sampleFromGaussianMixture(z_post_dens,component_mus,component_Cs,nSamp);
            % this might be very bad
            %samples = mvnrnd(repmat(mu_post',nSamp,1),C_post);
%             imagesc(samples)
%             pause
        end

        % transform values to firing rates
        rates = 10 * max(samples - 1,0).^(1.1);
        % calculate FR correlations
        C_r = cov(rates);
        
        % create spike counts
        windowLength = 20;
        nSpikeBin = floor(nSamp / windowLength);
        spikeCounts = zeros(nSpikeBin,size(C,1));
        for j=1:nSpikeBin
            actRate = sum(rates((j-1)*windowLength+1:j*windowLength,:),1);
            % rates are in Hz, a rate bin is 20 ms
            spikeCounts(j,:) = floor(actRate * (windowLength * 20 / 1000));
        end
        % calculate SC correlations
        C_s = cov(spikeCounts);
        
%         sum(all(C_post==0,2))
        sum(all(C_r==0,2))
               
        potCorr = upperTriangleValues(corrcov(C_post));
        rateCorr = upperTriangleValues(corrcov(C_r));     
        spikeCorr = upperTriangleValues(corrcov(C_s));     

        % plot everything
        subplot(3,nZ,i);
        hist(potCorr,linspace(-1,1,100));
        yl = ylim();
        if yl(2) > pmax
            pmax = yl(2);
        end
        subplot(3,nZ,nZ+i);
        hist(rateCorr,linspace(-1,1,100));
        yl = ylim();
        if yl(2) > rmax
            rmax = yl(2);
        end
        subplot(3,nZ,2*nZ+i);
        hist(spikeCorr,linspace(-1,1,100));
        yl = ylim();
        if yl(2) > smax
            smax = yl(2);
        end
    end
    for i=1:nZ
        subplot(3,nZ,i);
        title(sprintf('Z = %.2f',zs(i)),'FontSize',16);
        xlim([-1 1]);
        ylim([0 pmax*1.05]);
        if i==1
            ylabel('Potential','FontSize',16);
        end
        subplot(3,nZ,nZ+i);
        xlim([-1 1]);
        ylim([0 rmax*1.05]);
        if i==1
            ylabel('Rate','FontSize',16);
        end
        subplot(3,nZ,2*nZ+i);
        xlim([-1 1]);
        ylim([0 smax*1.05]);
        if i==1
            ylabel('Spike count','FontSize',16);
        end
    end
end