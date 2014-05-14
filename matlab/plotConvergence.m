function h = plotConvergence(ge,diffs,lsamples)
    if ischar(diffs)
        load(diffs);
    end
    diffs(:, ~any(diffs,1) ) = [];
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
    if lsamples == 0
        lrms = lumpedRMS(ge);
    else
        lrms = lumpedLike(ge,lsamples);
    end
    plot(xlim,[lrms(1) lrms(1)],'k--');
    plot(xlim,[lrms(2) lrms(2)],'r--');
    plot(xlim,[lrms(3) lrms(3)],'g--');
    %ylim([0,max(diffs(:,1))+0.2]);
    xlabel('Iterative EM step #');  
    if lsamples == 0
        ylabel('Mean squared error of covariance parameters')
    else
        ylabel('Unnormalised log-likelihood');
    end
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
    lrms(3) = 0;
end

function lrms = lumpedLike(ge,lsamples)
    lrms(1) = gestaltLogLikelihood(ge,lsamples);
    
    combined = ge.cc{1} + ge.cc{2};
    first = ge.cc{1};
    result = randomCovariances(2,ge.Dv);    
    
    result{1} = combined;
    ge.cc = result;
    lrms(2) = gestaltLogLikelihood(ge,lsamples);
    
    result{1} = first;
    ge.cc = result;
    lrms(3) = gestaltLogLikelihood(ge,lsamples);
end