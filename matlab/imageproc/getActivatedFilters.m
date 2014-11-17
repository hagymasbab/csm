function [activatedIndices,activatedCoeffs,handles] = getActivatedFilters(cc,A,coeffs,handles,plotFilters)
    actFilters = covariance2template(cc,A);
    activatedIndices = find(actFilters);
    activatedFilters = A(:,activatedIndices)';
    activatedCoeffs = coeffs(activatedIndices,1);
    
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
    
    if plotFilters
        if isempty(handles)
            handles(1) = figure('Units','normalized','OuterPosition',[0.05 0.1 0.4 0.8]);    
            handles(2) = figure('Units','normalized','OuterPosition',[0.55 0.1 0.4 0.8]);
        end
        figure(handles(1));
        viewImageSet(positiveFilters,'titles',positiveTitles);
        figure(handles(2));
        viewImageSet(negativeFilters,'titles',negativeTitles);
    end
end