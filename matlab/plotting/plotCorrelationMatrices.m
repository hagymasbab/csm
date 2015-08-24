function plotCorrelationMatrices()
    load('c_gsm-learned_248.mat')
    invprior = corrcov(inv(C));
    prior = corrcov(C);
    load('c_gsm-learnA-signal1000_248.mat')
    learnsignal = sc;
    load('c_gsm-synA-signal1000_248.mat')
    synsignal = sc;
    load('filters_gabor_256x248.mat')
    synfilter = corrcov(A'*A);
    invsynfilter = corrcov(stableInverse(A'*A));
    load('filters_gsm-learned_256x248.mat')
    learnfilter = corrcov(A'*A);
    invlearnfilter = corrcov(inv(A'*A));
    
    close all;
    
    subplot(2,2,1);
    scatter(upperTriangleValues(synfilter),upperTriangleValues(synsignal));
    xlabel('Gabor filter','FontSize',16);
    ylabel('Gabor signal','FontSize',16);
    
    subplot(2,2,2);
    scatter(upperTriangleValues(synfilter),upperTriangleValues(inv(synsignal)));
    xlabel('Gabor filter','FontSize',16);
    ylabel('Inverse Gabor signal','FontSize',16);
    
    subplot(2,2,3);
    scatter(upperTriangleValues(learnfilter),upperTriangleValues(learnsignal));
    xlabel('Learned filter','FontSize',16);
    ylabel('Learned signal','FontSize',16);
    
    subplot(2,2,4);
    scatter(upperTriangleValues(learnfilter),upperTriangleValues(inv(learnsignal)));
    xlabel('Learned filter','FontSize',16);
    ylabel('Inverse Learned signal','FontSize',16);
    
    
    

    figure;
    matrices = {synfilter,learnfilter,synsignal,learnsignal};
    names = {'Synthetic filter','Learned filter','Synthetic signal','Learned signal'};
    
    for i=1:4
        subplot(2,4,i);
        scatter(upperTriangleValues(prior),upperTriangleValues(matrices{i}));
        xlabel('Learned prior','FontSize',16);
        ylabel(names{i},'FontSize',16);
        subplot(2,4,4+i);
        scatter(upperTriangleValues(invprior),upperTriangleValues(matrices{i}));
        xlabel('Inverse learned prior','FontSize',16);
        ylabel(names{i},'FontSize',16);
    end    
end