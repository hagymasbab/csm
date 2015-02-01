function plotLikelihoodVariance(loadData,scientific,nTrials)
    %LS = [10 50 100 200 500 1000];
    %LS = [1000 500 200 50];
    LS = [10 100 1000];
    nls = length(LS);    
    %nTrials = 10;
    if ~loadData
        ge = gestaltCreate('temp','Dx',64,'k',2,'generateComponents',true,'generateData',true, ...
            'N',1,'filters','filters_OF_64.mat','nullComponent',false);
        lls = zeros(nls,nTrials);
        labels = cell(1,nls);
        for l = 1:nls     
            fprintf('NS %d ',LS(l));
            for t=1:nTrials
                printCounter(t,'maxVal',nTrials,'stringVal','trial');          
                lls(l,t) = gestaltLogLikelihood(ge,LS(l),ge.X,[],false,'shuffle',scientific);
            end
        end
        save('likelihood_varaiance.mat','lls');
    else
        load('likelihood_varaiance.mat');
    end
    
    means = zeros(nls,1);
    stds = zeros(nls,1);
    for l = 1:nls
        labels{l} = sprintf('%d',LS(l));
        act = lls(l,:)';
        act(isnan(act)) = [];
        means(l) = mean(act);
        stds(l) = std(act);
    end
    
    barwitherr(stds,-means);
    set(gca,'XTickLabel',labels,'FontSize',16);
    ylabel('Negative log-likelihood','FontSize',16);
    xlabel('Number of samples','FontSize',16);
end
        