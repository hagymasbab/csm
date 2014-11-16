function activatedCoeffs = plotActivatedFilters(cc,A,coeffs)
    actFilters = covariance2template(cc,A);
    cells = find(actFilters);
    activatedFilters = A(:,cells)';
    activatedCoeffs = coeffs(c,cells)';
    
    positiveFilters = activatedFilters(activatedCoeffs>0,:);
    negativeFilters = activatedFilters(activatedCoeffs<0,:);
    [positiveCoeffs,posperm] = sort(activatedCoeffs(activatedCoeffs>0,1),1,'descend');
    [negativeCoeffs,negperm] = sort(activatedCoeffs(activatedCoeffs<0,1),1,'ascend');        
    positiveFilters = positiveFilters(posperm,:);
    negativeFilters = negativeFilters(negperm,:);
    
    positiveTitles = {};
    negativeTitles = {};
    for i=1:length(positiveCoeffs)
        positiveTitles{i} = num2str(positiveCoeffs(i,1));
    end
    for i=1:length(negativeCoeffs)
        negativeTitles{i} = num2str(negativeCoeffs(i,1));
    end
    
    figure('Units','normalized','OuterPosition',[0.05 0.1 0.4 0.8]);
    viewImageSet(positiveFilters,'titles',positiveTitles);
    figure('Units','normalized','OuterPosition',[0.55 0.1 0.4 0.8]);
    viewImageSet(negativeFilters,'titles',negativeTitles);
end