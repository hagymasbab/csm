function gsmPosteriorCorrVsCont(A,C,sigma_x,zs,x_base)
    setrandseed(1);
    nZ = length(zs);
    if isempty(x_base)        
        iC = stableInverse(C);
        sATA = (1/sigma_x^2) * (A' * A);
        ymax = 0;
        for i=1:nZ
            actCov = stableInverse(iC + zs(i)^2 * sATA);
            actCorr = upperTriangleValues(corrcov(actCov));
            subplot(1,nZ,i);
            hist(actCorr,linspace(-1,1,100));
            yl = ylim();
            if yl(2) > ymax
                ymax = yl(2);
            end
        end
        for i=1:nZ
            subplot(1,nZ,i);
            xlim([-1 1]);
            ylim([0 ymax*1.05]);
        end
    else
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
        for i=1:nZ
            printCounter(i,'maxVal',nZ,'stringVal','Contrast');
            % construct x        
            x_act = zs(i) * x_base;
            % compute posterior
            [mu_post,C_post,~,~,z_post_dens,component_mus,component_Cs] = gsmPosteriorV(x_act,A,C,sigma_x,2,2,10);
            
            if any(isnan(z_post_dens))
                fprintf('NaN is Z dens\n');
                continue
            end
            % sample back from mixture
            %samples = sampleFromGaussianMixture(z_post_dens,component_mus,component_Cs,1000);
            % this might be very bad
            samples = mvnrnd(repmat(mu_post',1000000,1),C_post);
            imagesc(samples)
            pause
            
            % transform values to firing rates
            rates = 10 * max(samples - 1.9,0).^(1.1);

            % calculate FR correlations
            C_r = cov(rates);
            sum(all(C_post==0,2))
            sum(all(C_r==0,2))
            Cm = corrcov(C_r);
            rateCorr = upperTriangleValues(Cm);            
            potCorr = upperTriangleValues(corrcov(C_post));
            
            % plot everything
            subplot(2,nZ,i);
            hist(potCorr,linspace(-1,1,100));
            yl = ylim();
            if yl(2) > pmax
                pmax = yl(2);
            end
            subplot(2,nZ,nZ+i);
            hist(rateCorr,linspace(-1,1,100));
            yl = ylim();
            if yl(2) > rmax
                rmax = yl(2);
            end
        end
        for i=1:nZ
            subplot(2,nZ,i);
            xlim([-1 1]);
            ylim([0 pmax*1.05]);
            subplot(2,nZ,nZ+i);
            xlim([-1 1]);
            ylim([0 rmax*1.05]);
        end
    end
end