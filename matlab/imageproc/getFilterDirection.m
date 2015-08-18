function maxThetas = getFilterDirection(filters,thetaRes,loadStuff)

    nFilt = size(filters,2);
    imsize = sqrt(size(filters,1));
    lambdaRes = 5;
%     thetaRes = 180;
    phaseRes = 10;
    lambdaVals = linspace(2,imsize/2,lambdaRes);
    thetaVals = linspace(0,180,thetaRes);
    phaseVals = linspace(0,1,phaseRes);
        
    if loadStuff
        load('save_filterdir.mat');
    else
        maxThetas = zeros(nFilt,1);
        maxGrat = zeros(nFilt,imsize^2);
    end
    
    nPlot = min(10,nFilt);
    nRow = 2;
    plotIndices = [];
    for r = 1:nRow
        plotIndices = [plotIndices; chooseKfromN(nPlot,nFilt)];
    end
    toPlot = false(nFilt,1);
    toPlot(plotIndices(:)) = true;
    plotcount = 1;
    
    for i = 1:nFilt
        f = filters(:,i);
        if ~loadStuff
            printCounter(i,'maxVal',nFilt,'stringVal','Filter');            
            maxProd = -Inf;
            for l = 1:lambdaRes
                for p = 1:phaseRes
                    for t = 1:thetaRes
                        act_grat = grating(lambdaVals(l),thetaVals(t),phaseVals(p),imsize);
                        act_grat = act_grat(:);
                        actProd = f' * act_grat;
                        if actProd > maxProd
                            maxProd = actProd;
                            maxThetas(i) = thetaVals(t);
                            maxGrat(i,:) = act_grat';
                        end
                    end
                end
            end
        end
        if toPlot(i)
            actbundle = ceil(plotcount / nPlot);
            filterrow = (actbundle-1)* 2 + 1;
            gratingrow = actbundle * 2;
            actcol = plotcount - (actbundle-1) * nPlot;
            subplot(2*nRow,nPlot,(filterrow-1)*nPlot + actcol);
            viewImage(f,'useMax',true);
            title(sprintf('Filter no. %d',i));
            subplot(2*nRow,nPlot,(gratingrow-1)*nPlot + actcol);
            viewImage(maxGrat(i,:)');
            title(sprintf('Orientation %.2f',maxThetas(i)));
            plotcount = plotcount + 1;
        end
    end    
    
    if ~loadStuff
        save('bin/save_filterdir.mat','maxThetas','maxGrat');
    end
    
end