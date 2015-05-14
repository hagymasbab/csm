function correlationDispersion(cc,variance)
    k = length(cc);    
    colormap(bone);
    mean_cc = componentSum(1/k,cc);
    % choose the 20 largest values from all components
    nPairs = 10;
    valsPerComp = nPairs / k;
    maxcoords = zeros(nPairs,2);
    maxvals = zeros(nPairs,k);
    meanvals = zeros(nPairs,1);
    xlabels = {};
    
    def_colors = get(groot,'DefaultAxesColorOrder');
    blueish_color = def_colors(1,:);
    reddish_color = def_colors(2,:);
    
    plotrange = 5;
    excluded = [];
    for i=1:plotrange
        if variance
            cm = diag(diag(cc{i}));
        else
            cm = corrcov(cc{i});        
            cm = triu(abs(nodiag(cm)));
        end
        [~,xcoord,ycoord] = maxNElements(cm,valsPerComp,excluded);
        excluded = [excluded;xcoord ycoord];
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
            sc = scatter(idx*ones(k,1),maxvals(idx,:)',100,'filled');
            set(sc,'MarkerFaceColor',blueish_color);
            hold on;
            sc = scatter(idx,meanvals(idx,1),200,'filled','r');
            set(sc,'MarkerFaceColor',reddish_color);
        end                
    end
    xlim([0.5 plotrange+0.5]);
    %set(gca,'XTick',1:nPairs,'XTickLabel',xlabels);
    set(gca,'XTick',1,'XTickLabel',{''},'FontSize',20);
    %xlabel('Filter pairs','FontSize',24);
    %ylabel('Correlations in each component and mean','FontSize',16);
end

    