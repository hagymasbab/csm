function filterAndPriorCorrelations(learnedA,referenceA,C)
    nFilt = size(C,1);
    [maxProds,bestMatches] = max(referenceA'*learnedA);
    filterCorr = corrcov(learnedA'*learnedA);
    priorCorr = corrcov(C);
    %priorCorr = corrcov(C);
    shuffledPC = priorCorr(bestMatches,bestMatches);
    shuffledPC = zeros(nFilt);
    for i=1:nFilt
        for j=i+1:nFilt
            shuffledPC(i,j) = priorCorr(bestMatches(i),bestMatches(j));
        end
    end
    
%     subplot(1,2,1)
    scatter(upperTriangleValues(nodiag(priorCorr)),upperTriangleValues(nodiag(shuffledPC)));
    corrcoef(upperTriangleValues(nodiag(priorCorr)),upperTriangleValues(nodiag(shuffledPC)))
    length(upperTriangleValues(nodiag(priorCorr)))
    xlabel('Prior correlations of learned filters','FontSize',16);
    ylabel('Prior correlations of matched filters','FontSize',16);
%     
%     subplot(1,2,2)
%     scatter(upperTriangleValues(nodiag(filterCorr)),upperTriangleValues(nodiag(priorCorr)));
%     xlabel('Learned filter correlations','FontSize',16);
%     ylabel('Inverse prior correlations','FontSize',16);
end
        