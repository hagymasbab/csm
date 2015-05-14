function contrastDistributions(N,loadSamples,randseed)
    
    close all;
    setrandseed(randseed);
    k = 2;
    nSamples = 50;
    burnin = 20;
    
    Dx = [256 576];
    Z = zeros(length(Dx),N);
    xmin = Inf;
    xmax = -Inf;
    ymax = -Inf;
    for x = 1:length(Dx)
        if loadSamples
            load(sprintf('contr_dist_samples_%d.mat',x));
        else
            load(sprintf('patches_vanhateren_%d.mat',Dx(x)));
            stim_idx = chooseKfromN(N,size(patchDB,2));
            stimuli = cell(1,N);
            for i=1:N
                stimuli{i} = patchDB(:,stim_idx(i));
            end

            timings = (nSamples+burnin) * ones(1,N); 
            ge = gestaltCreate('temp','Dx',Dx(x),'k',k,'B',1,'N',1,'filters','OF','obsVar',1,'z_shape',1, ...
                'nullComponent',false,'generateComponents',true,'generateData',false);

            [~,~,zsamp] = gestaltScheduling(stimuli,timings,{ge},1,ge.obsVar,true,'gibbs');
            save(sprintf('contr_dist_samples_%d.mat',x),'zsamp');
        end
        zdata = squeeze(zsamp(1,1,:));
        for i=1:N
            start_idx = (i - 1) * (nSamples+burnin) + burnin + 1;
            end_idx = start_idx + nSamples - 1;
            Z(x,i) = mean(zdata(start_idx:end_idx));
        end
        subplot(length(Dx),1,x);
        hist(Z(x,:));
        title(sprintf('Estimated contrasts for %d images, D_x = %d',N,Dx(x)),'FontSize',16);
        xl = xlim();
        yl = ylim();
        if xl(1)<xmin, xmin = xl(1); end
        if xl(2)>xmax, xmax = xl(2); end
        if yl(2)>ymax, ymax = yl(2); end        
    end
    for i = 1:length(Dx)
        subplot(length(Dx),1,i);
        xlim([xmin xmax]);
        ylim([0 ymax]);
        set(gca, 'FontSize',16);
    end
end
        