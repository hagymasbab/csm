function gsmOrientSelEM(iterfile,loadStuff,posterior,thetaRes,FRnonlin)
    setrandseed(1);
    load(iterfile);
    nStep = length(A_iter);
    nFilt = size(A_iter{1},2);
    imsize = sqrt(size(A_iter{1},1));
    if loadStuff
        load('save_emorient.mat');
    else        
        tuning_curves = zeros(nStep,nFilt,thetaRes);
        respVar = zeros(nStep,nFilt);
        prefThetas = zeros(nStep,nFilt);
        for i=1:nStep        
            [prefTheta,tuningCurves,prefPhase,prefLambda,prefVars] = gsmOrientationSelectivity(A_iter{i},C,sigma_iter(i),thetaRes,false,'grating',~posterior,1,false);
            prefThetas(i,:) = prefTheta;
            %pause
            tuning_curves(i,:,:) = tuningCurves;
            if posterior
                respVar(i,:) = prefVars;
            end
            for j = 1:nFilt                                
                if ~posterior
%                     printCounter(j,'StringVal','Preferred posterior','maxVal',nFilt);
%                     prefGrat = grating(prefLambda(j),prefTheta(j),prefPhase(j),imsize);
%                     [~,C_post,~,~] = gsmPosteriorV(prefGrat(:),A_iter{i},C,sigma_iter(i),2,2,10);
%                     respVar(i,j) = C_post(j,j);
                end
            end
        end
        save('bin/save_emorient.mat','respVar','tuning_curves','prefThetas');
    end
    
    OSIs = zeros(nStep,nFilt);
    for i=1:nStep        
        for j = 1:nFilt                
            OSIs(i,j) = orientationSelectivityIndex(tuning_curves(i,j,:),FRnonlin);
        end
        hist(prefThetas(i,:),linspace(0,180,100));
        title(sprintf('EM step %d',em_steps(i)))
        pause
    end
    
    close all;
    errorbar(em_steps,mean(OSIs,2),std(OSIs,0,2),'LineWidth',2);
    xlim([0 em_steps(nStep)+20]);
    xlabel('EM step #','FontSize',16);
    ylabel('OSI mean and st.d.','FontSize',16);
%     figure;
%     errorbar(em_steps,mean(respVar,2),std(respVar,0,2),'LineWidth',2);
%     xlim([0 em_steps(nStep)+20]);
%     xlabel('EM step #','FontSize',16);
%     ylabel('Response variance to pref. stim., mean and st.d.','FontSize',16);
end
        