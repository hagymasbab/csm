function learnedLikelihoodVariance(ge,cc,N_test,nTrial,samplenums,loadData,plotStuff)
    load('patches_vanhateren_576.mat'); % TODO gerealise
    N_all = size(patchDB,2);
    X_test = patchDB(:,chooseKfromN(N_test,N_all))';
    ge.cc = cc;
    if loadData
        load('bin/save_likevar.mat');
    else
        ll = zeros(nTrial,length(samplenums));
        for s = 1:length(samplenums)
            parfor t = 1:nTrial
                ll(t,s) = gestaltLogLikelihood(ge,samplenums(s),X_test,'scientific',true); 
            end
        end
        save('bin/save_likevar.mat','ll');
    end
    
    if plotStuff
        barwitherr(std(ll)',mean(ll)');
        xlabels={};for i=1:length(samplenums);xlabels{end+1}=sprintf('%d',samplenums(i));end;
        set(gca,'XTickLabel',xlabels)
    end
end