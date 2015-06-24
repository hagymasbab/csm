function values = plotFunctionProjections(func,resolution,mincoord,maxcoord)
    % plot 3 contours of 2D projections of a 3D function
    xgrid = linspace(mincoord,maxcoord,resolution);
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
    proj1 = sum(proj12,2);
    proj2 = sum(proj12,1);
    proj3 = sum(proj23,1);
    
    clines = 30;
    
    clf;
    subplot(2,3,1);    
    plot(xgrid,proj1);    
    subplot(2,3,2);
    plot(xgrid,proj2);
    subplot(2,3,3);
    plot(xgrid,proj3);
    
    subplot(2,3,4);
    contourf(xgrid,xgrid,proj12,clines);
    %surf(xgrid,xgrid,proj12);
    subplot(2,3,5);
    contourf(xgrid,xgrid,proj13,clines);
    subplot(2,3,6);
    contourf(xgrid,xgrid,proj23,clines);
    
    %figureHandle = inputFigureHandle;
    
end