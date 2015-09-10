function gsmPosteriorCorrVsCont(A,C,sigma_x,zs,x_base,rateTh)
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
    
    nRow = 4;
    maxvals = zeros(nRow,1);
    nSamp = 20000;
    for i=1:nZ
        printCounter(i,'maxVal',nZ,'stringVal','Contrast');

        if isempty(x_base)        
            C_post = stableInverse(iC + zs(i)^2 * sATA);
            samples = mvnrnd(zeros(nSamp,size(A,2)),C_post);
            mean(abs(samples(:)))
        else
            % construct x        
            x_act = zs(i) * x_base;
            % compute posterior
            [mu_post,C_post,~,~,z_post_dens,component_mus,component_Cs] = gsmPosteriorV(x_act,A,C,sigma_x,2,2,10);
            mean(mu_post)
            mean(abs(mu_post))

            if any(isnan(z_post_dens))
                fprintf('NaN is Z dens\n');
                continue
            end
            % sample back from mixture
            samples = sampleFromGaussianMixture(z_post_dens,component_mus,component_Cs,nSamp);
            mean(abs(samples(:)))
            % this might be very bad
            %samples = mvnrnd(repmat(mu_post',nSamp,1),C_post);
%             imagesc(samples)
%             pause
        end
        
        C_samp = cov(samples);
        
        % transform values to firing rates
        rates = 10 * max(samples - rateTh,0).^(1.1);
        %rates = 10 * max(abs(samples) - rateTh,0).^(1.1);
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
        potCorr(isnan(potCorr)) = 0;
        sampCorr = upperTriangleValues(corrcov(C_samp));     
        sampCorr(isnan(sampCorr)) = 0;
        rateCorr = upperTriangleValues(corrcov(C_r));     
        rateCorr(isnan(rateCorr)) = 0;
        spikeCorr = upperTriangleValues(corrcov(C_s));     
        spikeCorr(isnan(spikeCorr)) = 0;                
        allcorrs = [potCorr sampCorr rateCorr spikeCorr];   
        names = {'Potential','Sample','Rate','Spike count'};
                
        % plot everything
        for r = 1:nRow
            subplot(nRow,nZ,(r-1)*nZ + i);
            hist(allcorrs(:,r),linspace(-1,1,100));
            yl = ylim();
            if yl(2) > maxvals(r)
                maxvals(r) = yl(2);
            end
        end
%         subplot(nRow,nZ,nZ+i);
%         hist(rateCorr,linspace(-1,1,100));
%         yl = ylim();
%         if yl(2) > rmax
%             rmax = yl(2);
%         end
%         subplot(3,nZ,2*nZ+i);
%         hist(spikeCorr,linspace(-1,1,100));
%         yl = ylim();
%         if yl(2) > smax
%             smax = yl(2);
%         end
    end
    for i=1:nZ
        for r = 1:nRow
            subplot(nRow,nZ,(r-1)*nZ + i);
            xlim([-1 1]);
            ylim([0 maxvals(r)*1.05]);
            if i==1
                ylabel(names{r},'FontSize',16);
            end
        end
        
%         subplot(3,nZ,i);
%         title(sprintf('Z = %.2f',zs(i)),'FontSize',16);
        
%         xlim([-1 1]);
%         ylim([0 pmax*1.05]);
%         if i==1
%             ylabel('Potential','FontSize',16);
%         end
%         subplot(3,nZ,nZ+i);
%         xlim([-1 1]);
%         ylim([0 rmax*1.05]);
%         if i==1
%             ylabel('Rate','FontSize',16);
%         end
%         subplot(3,nZ,2*nZ+i);
%         xlim([-1 1]);
%         ylim([0 smax*1.05]);
%         if i==1
%             ylabel('Spike count','FontSize',16);
%         end
    end
end