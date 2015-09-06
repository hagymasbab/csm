function [maxtun,responses,maxPhase,maxLambda,prefRespVars] = gsmOrientationSelectivity(A,C,sigma_x,thetaRes,loadStuff,match,scalarprod,contrast,plotStuff)

    setrandseed(1);
    nFilt = size(A,2);
    imsize = sqrt(size(A,1));
    lambdaRes = 3;
    phaseRes = 10;
    spaceRes = 4;
    lambdaVals = linspace(2,imsize/2,lambdaRes);
    thetaVals = linspace(0,180,thetaRes);
    phaseVals = linspace(0,1,phaseRes);
    xVals = round(linspace(1,imsize,spaceRes));
    yVals = round(linspace(1,imsize,spaceRes));
        
    if loadStuff
        load('save_gsmorient.mat');
    else
        responses = zeros(nFilt,thetaRes);  
        respVars = zeros(nFilt,thetaRes);  
        respLambda = zeros(nFilt,thetaRes);
        respPhase = zeros(nFilt,thetaRes);
    end
    
    if plotStuff
        nPlot = min(10,nFilt);
        nRow = 2;
        plotIndices = [];
        for r = 1:nRow
            plotIndices = [plotIndices; chooseKfromN(nPlot,nFilt)];
        end
        toPlot = false(nFilt,1);
        toPlot(plotIndices(:)) = true;
        plotcount = 1;
    end
    
%     z = 1;
%     gsmPostMeanTransform = (z / sigma_x^2) * stableInverse(stableInverse(C) + (z^2 / sigma_x^2) * (A'*A)) * A';
        
    if ~loadStuff         
        for t = 1:thetaRes
            printCounter(t,'maxVal',thetaRes,'stringVal','Orientation');    
            theta_responses = zeros(lambdaRes,nFilt);
            theta_variances = zeros(lambdaRes,nFilt);
            phase_indices = zeros(lambdaRes,nFilt);
            for l = 1:lambdaRes                
                if strcmp(match,'grating')
                    lambda_responses = zeros(phaseRes,nFilt);
                    lambda_variances = zeros(phaseRes,nFilt);
                    for p = 1:phaseRes                    
                        act_grat = contrast * grating(lambdaVals(l),thetaVals(t),phaseVals(p),imsize);                            
                        act_stim = act_grat(:);
                        %act_stim = act_stim ./ std(act_stim);
                        %act_resp = gsmPostMeanTransform * act_stim;
                        if scalarprod
                            act_resp = A' * act_stim;
                            act_var = zeros(nFilt,1);
                        else
                            %act_resp = gsmPostMeanTransform * act_stim;
                            [act_resp,act_C_post] = gsmPosteriorV(act_stim,A,C,sigma_x,2,2,10);
                            act_var = diag(act_C_post);
                        end
                        %responses(:,t) = responses(:,t) + act_resp;
                        lambda_responses(p,:) = act_resp;   
                        lambda_variances(p,:) = act_var;  
                    end                    
                elseif strcmp(match,'gabor')
                    lambda_responses = zeros(spaceRes^2,nFilt);
                    % tODO lambda variances
                    for x = 1:spaceRes
                        for y = 1:spaceRes
                            act_gabor = contrast * gaborfilter(lambdaVals(l),thetaVals(t),imsize,xVals(x),yVals(y));
                            act_stim = act_gabor(:);
                            act_stim = act_stim ./ std(act_stim);
                            %act_resp = gsmPostMeanTransform * act_stim;
                            %responses(:,t) = responses(:,t) + act_resp;                            
                            if scalarprod
                                act_resp = A' * act_stim;
                            else
                                act_resp = gsmPostMeanTransform * act_stim;
                            end
                            lambda_responses((x-1)*spaceRes+y,:) = act_resp;
                        end
                    end
                else
                    error('no valid match');
                end
                % this only makes sense if we used gratings
                [theta_resp,phase_idx] = max(lambda_responses);  
                theta_responses(l,:) = theta_resp';
                theta_variances(l,:) = lambda_variances(phase_idx);
                phase_indices(l,:) = phase_idx';
            end
            [resp,lambda_idx] = max(theta_responses);
            responses(:,t) = resp';
            respVars(:,t) = theta_variances(lambda_idx);
            phase_idx = phase_indices(lambda_idx');
            actPhase = phaseVals(phase_idx);
            actLambda = lambdaVals(lambda_idx');
            respPhase(:,t) = actPhase;
            respLambda(:,t) = actLambda;
        end
    end
    
    [~,maxidx] = max(responses');
    maxtun = thetaVals(maxidx);
    maxPhase = respPhase(maxidx);
    maxLambda = respLambda(maxidx);
    trespVars = respVars';
    prefRespVars = trespVars(maxidx);
    
    if plotStuff
        close all;
        priorcorr_vs_orientdiff = [];
        priorCorr = corrcov(C);
        for i = 1:nFilt
            for j=i+1:nFilt
                priorcorr_vs_orientdiff = [priorcorr_vs_orientdiff; abs(maxtun(i)-maxtun(j)) priorCorr(i,j)];
            end
        end
        scatter(priorcorr_vs_orientdiff(:,1),abs(priorcorr_vs_orientdiff(:,2)));
        xlabel('Preferred orientation difference','FontSize',16);
        ylabel('Prior correlation magnitude','FontSize',16);
%         hold on
%         correlationPlot(priorcorr_vs_orientdiff(:,1),abs(priorcorr_vs_orientdiff(:,2)));
%         hold off
        figure;
    
        for i = 1:nFilt
            f = A(:,i);
            if toPlot(i)
                rownum = 2;
                actbundle = ceil(plotcount / nPlot);
                filterrow = (actbundle-1)* rownum + 1;
                tuningrow = actbundle * rownum;
                actcol = plotcount - (actbundle-1) * nPlot;

                subplot(rownum*nRow,nPlot,(filterrow-1)*nPlot + actcol);
                viewImage(f,'useMax',true);
                title(sprintf('Filter no. %d',i));            

                subplot(rownum*nRow,nPlot,(tuningrow-1)*nPlot + actcol);
                plot(thetaVals,responses(i,:))
                [~,maxidx] = max(responses(i,:));
                title(sprintf('Max tuning %.2f',maxtun(i)));

                plotcount = plotcount + 1;
            end
        end    
    end
    
    if ~loadStuff
        save('bin/save_gsmorient.mat','responses','respLambda','respPhase','respVars');
    end
    
end

function img = gaborfilter(lambda,theta,imSize,xc,yc) 
    thetaRad = (theta / 360) * 2*pi;
    lambda = lambda - 1;
    px = xc/imSize;
    py = yc/imSize;
    [~,~,img] = gabor('theta',thetaRad,'lambda',lambda,'width',imSize,'height',imSize,'px',px,'py',py,'Sigma',lambda/2);             
end
