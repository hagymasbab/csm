function [values,figureHandle] = plotFunctionProjections(func,stepsize,inputFigureHandle)
    % plot 3 contours of 2D projections of a 3D function
    xgrid = 0:stepsize:1;
    dimsize = size(xgrid,2);
    values = zeros(dimsize,dimsize,dimsize);
    for i=1:dimsize
        for j=1:dimsize
            for k=1:dimsize
                values(i,j,k) = func(xgrid(i),xgrid(j),xgrid(k));
            end
        end
    end
    
    proj12 = squeeze(sum(values,3))';
    proj13 = squeeze(sum(values,2))';
    proj23 = squeeze(sum(values,1))';
    
    clines = 10;
    
    if isnan(inputFigureHandle)
        figureHandle = figure();
    else
        figure(inputFigureHandle);
        figureHandle = inputFigureHandle;
    end
    clf;
    subplot(1,3,1);
    contour(xgrid,xgrid,proj12,clines);
    subplot(1,3,2);
    contour(xgrid,xgrid,proj13,clines);
    subplot(1,3,3);
    contour(xgrid,xgrid,proj23,clines);
    
end