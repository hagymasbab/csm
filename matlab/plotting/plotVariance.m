function vseries = plotVariance(nz,randseed)
    if strcmp(randseed,'last')
        load lastrandseed;
    end
    s = RandStream('mt19937ar','Seed',randseed);
    RandStream.setGlobalStream(s);
    randseed = s.Seed;
    save('lastrandseed.mat','randseed');

    Dv = 100;
    nSamp = 200;
    burnin = 100;
    timings = ones(1,nz) * nSamp;
    [vdata,zdata,vsamp] = gestaltPriorCorrelations(1,timings,'','transient','contrast','priors','gibbs',Dv);
    zdata = squeeze(zdata(1,1,1,:));
    vdata = squeeze(vdata(1,1,1,:));
    
    vseries = reshape(vsamp(1,1,:,1,:),nz*nSamp,Dv);
    %vseries = vseries(:,burnin+1:end,:);
    %variances = squeeze(std(vseries,0,2));
    variances = [];
    zlabel = {};
    for i=1:nz
        actvar = [];
        for j = 1:Dv
            actvar = [actvar std(vseries((i-1)*nSamp+1:i*nSamp,j))];
        end
        actz = mean(zdata((i-1)*nSamp+1:i*nSamp));
        zlabel{end+1} = sprintf('%.2f',actz);
        variances = [variances; actvar];
    end
    mean_vars = mean(variances,2);
    std_vars = std(variances,0,2);
    figure
    barwitherr(std_vars,mean_vars);
    xlabel('estimated contrast');
    ylabel('V1 response variance');
    set(gca,'XTickLabel',zlabel);
    set(gca,'FontSize',24);
    xlim([0.5 nz+0.5]);
    yl = ylim();
    ylim([0 yl(2)]);
    
    kurt = [];
    for i=1:nz*nSamp
        kurt = [kurt; kurtosis(vseries(i,:)') - 3];
    end
    figure
    scatter(zdata,kurt);
    xlabel('contrast');
    ylabel('V1 excess kurtosis');
    set(gca,'FontSize',24);
    xlim([0 2]);
    yl = ylim();
    ylim([-1 yl(2)]);
    
end
    