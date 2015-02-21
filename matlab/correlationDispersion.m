function correlationDispersion(cc)
    k = length(cc);
    mean_cc = componentSum(1/k,cc);
    % choose the 20 largest values from all components
    nPairs = 20;
    valsPerComp = nPairs / k;
    maxcoords = zeros(nPairs,2);
    maxvals = zeros(nPairs,k);
    meanvals = zeros(nPairs,1);
    xlabels = {};
    for i=1:k
        cm = corrcov(cc{i});
        cm = triu(abs(nodiag(cm)));
        [~,xcoord,ycoord] = maxNElements(cm,valsPerComp);
        start_idx = (i-1)*valsPerComp + 1;
        end_idx = start_idx + valsPerComp - 1;
        maxcoords(start_idx:end_idx,:) = [xcoord ycoord];
        
        for h=1:valsPerComp            
            idx = (i-1)*valsPerComp+h;
            xlabels{end+1} = sprintf('%d %d',maxcoords(idx,1),maxcoords(idx,2));
            for j=1:k                
                maxvals(idx,j) = cc{j}(maxcoords(idx,1),maxcoords(idx,2));
            end
            meanvals(idx,1) = mean(maxvals(idx,:),2);
            scatter(idx*ones(k,1),maxvals(idx,:)',[],'b')
            hold on;
            scatter(idx,meanvals(idx,1),'filled','r')
        end                
    end
    xlim([0 nPairs+1]);
    set(gca,'XTick',1:nPairs,'XTickLabel',xlabels);
    xlabel('Filter pair','FontSize',16);
    ylabel('Correlations in each component and mean','FontSize',16);
end

    