function maxThetas = getFilterDirection(filters,thetaRes,loadStuff,match)

    nFilt = size(filters,2);
    imsize = sqrt(size(filters,1));
    lambdaRes = 3;
%     thetaRes = 180;
    phaseRes = 10;
    spaceRes = 4;
    lambdaVals = linspace(2,imsize/2,lambdaRes);
    thetaVals = linspace(0,180,thetaRes);
    phaseVals = linspace(0,1,phaseRes);
    xVals = round(linspace(1,imsize,spaceRes));
    yVals = round(linspace(1,imsize,spaceRes));
        
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
    
    %[matched_gabors,gabor_perm] = matchGabors(filters,'filters_gabor_256x248.mat',1,false);
    
    for i = 1:nFilt
        f = filters(:,i);
        %gf = matched_gabors(:,i);
        if ~loadStuff
            printCounter(i,'maxVal',nFilt,'stringVal','Filter');            
            maxProd = -Inf;
            for l = 1:lambdaRes
                for t = 1:thetaRes
                    if strcmp(match,'grating')
                        for p = 1:phaseRes                    
                            act_grat = grating(lambdaVals(l),thetaVals(t),phaseVals(p),imsize);
                            act_grat = act_grat(:);
                            actProd = f' * act_grat;
                            if actProd > maxProd
                                maxProd = actProd;
                                maxThetas(i) = thetaVals(t);
                                maxGrat(i,:) = act_grat';
                            end
                        end
                    elseif strcmp(match,'gabor')
                        for x = 1:spaceRes
                            for y = 1:spaceRes
                                act_gabor = gaborfilter(lambdaVals(l),thetaVals(t),imsize,xVals(x),yVals(y));
                                act_gabor = act_gabor(:);
                                actProd = f' * act_gabor;
                                if actProd > maxProd
                                    maxProd = actProd;
                                    maxThetas(i) = thetaVals(t);
                                    maxGrat(i,:) = act_gabor';
                                end
                            end
                        end
                    else
                        error('no valid match');
                    end
                end
            end
        end
        if toPlot(i)
            rownum = 2;
            actbundle = ceil(plotcount / nPlot);
            filterrow = (actbundle-1)* rownum + 1;
            %gaborrow = (actbundle-1)* rownum + 2;
            gratingrow = actbundle * rownum;
            actcol = plotcount - (actbundle-1) * nPlot;
            
            subplot(rownum*nRow,nPlot,(filterrow-1)*nPlot + actcol);
            viewImage(f,'useMax',true);
            title(sprintf('Filter no. %d',i));
            
%             subplot(rownum*nRow,nPlot,(gaborrow-1)*nPlot + actcol);
%             viewImage(gf,'useMax',true);
%             title(sprintf('Gabor no. %d',gabor_perm(i)));
            
            subplot(rownum*nRow,nPlot,(gratingrow-1)*nPlot + actcol);
            viewImage(maxGrat(i,:)');
            title(sprintf('Orientation %.2f',maxThetas(i)));
            
            plotcount = plotcount + 1;
        end
    end    
    
    if ~loadStuff
        save('bin/save_filterdir.mat','maxThetas','maxGrat');
    end
    
end

function img = gaborfilter(lambda,theta,imSize,xc,yc) 
    thetaRad = (theta / 360) * 2*pi;
    lambda = lambda - 1;
    px = xc/imSize;
    py = yc/imSize;
    [~,~,img] = gabor('theta',thetaRad,'lambda',lambda,'width',imSize,'height',imSize,'px',px,'py',py,'Sigma',lambda/2);             
end