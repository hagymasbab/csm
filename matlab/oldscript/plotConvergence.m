function h = plotConvergence(ge,diffs,lsamples)
    if ischar(diffs)
        load(diffs);
    end
    %diffs(:, ~any(diffs,1) ) = [];
    %diffs(diffs == 0) = [];
    for i = 1:size(diffs,1)
       for j = 2:size(diffs,2) 
           if diffs(i,j) == 0
               diffs(i,j) = diffs(i,j-1);
           end
       end
    end
    legends = cell(1,size(diffs,1));
    %legends{1} = 'mean & 1';
    for i=1:size(diffs,1)
        legends{i} = int2str(i);
    end
    h=figure();
    ha = axes();
    errorbar(0:size(diffs,2)-1,mean(diffs,1),std(diffs,1),'Parent',ha);    
    set( findobj(gca,'type','line'), 'LineWidth', 3);
    hold on;
    p1 = plot(0:size(diffs,2)-1,diffs,'Parent',ha);   
    leg1 = legend(p1,legends,'Location','NorthEastOutside');
    
    if lsamples == 0
        lrms = lumpedRMS(ge);
    else
        %lrms = lumpedLike(ge,lsamples);
        lrms = zeros(1,3);
    end
    %plot(xlim,[lrms(1) lrms(1)],'k--');
    %plot(xlim,[lrms(2) lrms(2)],'r--');
    %plot(xlim,[lrms(3) lrms(3)],'g--');
    p2 = plot(repmat(xlim,3,1)',[lrms' lrms']','--');    
    xlabel('Iterative EM step #');  
    if lsamples == 0
        ylim([0,max(diffs(:,1))+0.2]);
        ylabel('Mean squared error of covariance parameters')
    else
        ylabel('Unnormalised log-likelihood');
    end
    
    leg2 = copyobj(leg1,h);
    %delete( get(ah,'Children') );
    %set(hAx(2), 'Color','none', 'XTick',[],'YAxisLocation','right', 'Box','off')
    legend(p2,{'true','combined','single'},'Location','SouthEastOutside');
end

function lrms = lumpedRMS(ge)
    % only works for 2 components, covariance formulation and lumping into
    % the first
    trueSamples = reshape([ge.G reshape(ge.V,ge.N,ge.Dv*ge.B)],ge.N,1,ge.k+ge.Dv*ge.B);
    gVV = weightedSampleCovariance(trueSamples,ge.k,ge.B);
    %ge.cc = gVV;
    lrms(1) = covcompRootMeanSquare(gVV,ge.cc,[1 2]);
    combined = gVV{1} + gVV{2};
    result = randomCovariances(2,ge.Dv);
    result{1} = combined;
    lrms(2) = covcompRootMeanSquare(result,ge.cc,[1 2]);
    result{1} = gVV{1};
    lrms(3) = covcompRootMeanSquare(result,ge.cc,[1 2]);
    
end

function lrms = lumpedLike(ge,lsamples)
    trueSamples = reshape([ge.G reshape(ge.V,ge.N,ge.Dv*ge.B)],ge.N,1,ge.k+ge.Dv*ge.B);
    gVV = weightedSampleCovariance(trueSamples,ge.k,ge.B);
    
    ge.cc = gVV;    
    lrms(1) = gestaltLogLikelihood(ge,lsamples,0);

    rands = randomCovariances(2,ge.Dv);    
    
    ge.cc{1} = gVV{1} + gVV{2};
    ge.cc{2} = rands{2};
    lrms(2) = gestaltLogLikelihood(ge,lsamples,0);   
    
    ge.cc{1} = gVV{1};
    lrms(3) = gestaltLogLikelihood(ge,lsamples,0);
end