function signalAndPriorCorrelations(A,C,perm)
    numFilt = size(A,2);
    priorCorr = corrcov(C);
%     %priorCorr = log(abs(priorCorr));
%     priorCorr = log(priorCorr.^2);
%     %priorCorr = 1 ./ priorCorr;
    filterCorr = A'*A;
%     %filterCorr = log(abs(filterCorr));
%     filterCorr = log(filterCorr.^2);
    
    corrPairs = [];
    for i=1:numFilt
        for j=i+1:numFilt
            if perm(i) ~= perm(j)
                corrPairs = [corrPairs; filterCorr(i,j) priorCorr(perm(i),perm(j))];
            end
        end
    end
    
    scatter(corrPairs(:,1),corrPairs(:,2));
    cm = corr(corrPairs);
    hold on;
    correlationPlot(corrPairs(:,1),corrPairs(:,2),cm(1,2));
    hold off;
end