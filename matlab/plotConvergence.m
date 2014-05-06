function h = plotConvergence(ge,diffs)
    if ischar(diffs)
        load(diffs);
    end
    legends = cell(1,size(diffs,1)+1);
    legends{1} = 'mean';
    for i=1:size(diffs,1)
        legends{i+1} = int2str(i);
    end
    h=figure();
    errorbar(0:size(diffs,2)-1,mean(diffs,1),std(diffs,1));
    set( findobj(gca,'type','line'), 'LineWidth', 3);
    hold on;
    plot(0:size(diffs,2)-1,diffs);   
    lrms = lumpedRMS(ge);
    plot(xlim,[lrms(1) lrms(1)],'k--');
    plot(xlim,[lrms(2) lrms(2)],'k--');
    ylim([0,max(diffs(:,1))+0.2]);
    xlabel('Iterative EM step #');    
    ylabel('Mean squared error of covariance parameters')
    legend(legends,'Location','NorthEastOutside');
end

function lrms = lumpedRMS(ge)
    % only works for 2 components, covariance formulation and lumping into
    % the first
    combined = ge.cc{1} + ge.cc{2};
    result = randomCovariances(2,ge.Dv);
    result{1} = combined;
    lrms(1) = covcompRootMeanSquare(result,ge.cc,[1 2]);
    result{1} = ge.cc{1};
    lrms(2) = covcompRootMeanSquare(result,ge.cc,[1 2]);
end