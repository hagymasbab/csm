function plotLikelihoodVariance(loadData,nTrials,N,randseed)
    setrandseed(randseed);
    %LS = [10 50 100 200 500 1000];
    %LS = [1000 500 200 50];
    LS = [10 20 50 100 200];
    nls = length(LS);    
    %nTrials = 10;
    if ~loadData
        ge = gestaltCreate('temp','Dx',64,'k',2,'generateComponents',true,'generateData',true, ...
            'N',N,'filters','OF','nullComponent',false,'componentShape','oriented-gabors');
        %lls_sci = zeros(nls,nTrials);
        lls_nosci = zeros(nls,nTrials);
        labels = cell(1,nls);
        for l = 1:nls     
            fprintf('NS %d ',LS(l));
            for t=1:nTrials
                printCounter(t,'maxVal',nTrials,'stringVal','trial');          
                %act_ll_sci = gestaltLogLikelihood(ge,LS(l),ge.X,'scientific',true);
                act_ll_nosci = gestaltLogLikelihood(ge,LS(l),ge.X);
                %lls_sci(l,t) = act_ll_sci;
                lls_nosci(l,t) = act_ll_nosci;
                %pause
            end
        end
        save('likelihood_varaiance.mat','lls_sci','lls_nosci');
    else
        load('likelihood_varaiance.mat');
    end
    
    %means_sci = zeros(nls,1);
    %stds_sci = zeros(nls,1);
    means_nosci = zeros(nls,1);
    stds_nosci = zeros(nls,1);
    
    for l = 1:nls
        labels{l} = sprintf('%d',LS(l));
        
        %act = lls_sci(l,:)';
        %act(isnan(act)) = [];
        %means_sci(l) = mean(act);
        %stds_sci(l) = std(act);
        
        act = lls_nosci(l,:)';
        act(isnan(act)) = [];
        means_nosci(l) = mean(act);
        stds_nosci(l) = std(act);
    end
    
    %allstds = [stds_sci stds_nosci];
    allstds = stds_nosci;
    %allmeans = [-means_sci -means_nosci];
    allmeans = -means_nosci;
    
    %barwitherr([stds_sci stds_nosci],[-means_sci -means_nosci]);
    barwitherr(allstds,allmeans);
    colormap summer;
    ymin = min(allmeans(:)) - max(allstds(:)) - 3;
    ymax = max(allmeans(:)) + max(allstds(:)) + 3;
    ylim([ymin ymax]);
    set(gca,'XTickLabel',labels,'FontSize',16);
    ylabel('Negative log-likelihood','FontSize',16);
    xlabel('Number of samples','FontSize',16);
end
        