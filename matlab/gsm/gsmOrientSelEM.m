function gsmOrientSelEM(iterfile,loadStuff,posterior,thetaRes,FRnonlin,plotPhases)
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
        zmeans = zeros(nStep,nFilt);
        zvars = zeros(nStep,nFilt);
        for i=1:nStep       
            [prefTheta,tuningCurves,prefPhase,prefLambda,prefVars,zmean,zvar] = gsmOrientationSelectivity(A_iter{i},C,sigma_iter(i),thetaRes,false,'grating',~posterior,1,plotPhases); 
            zmeans(i,:) = zmean;
            zvars(i,:) = zvar;
            prefThetas(i,:) = prefTheta;
            if plotPhases
                pause
            end
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
        save('bin/save_emorient.mat','respVar','tuning_curves','prefThetas','zmeans','zvars');
    end
    
    close all;
    OSIs = zeros(nStep,nFilt);
    postVars = zeros(nStep,nFilt);
    for i=1:nStep        
        load('testX.mat')
        %x = randn(size(A_iter{i},1),1);
        [act_resp,act_C_post,z_mean,z_var] = gsmPosteriorV(x',A_iter{i},C,sigma_iter(i),2,2,20);
        postVars(i,:) = diag(act_C_post)';
        
        ATA = A_iter{i}' * A_iter{i};
        iATA = stableInverse(ATA);
        
        subplot(4,nStep,i);
        hist(upperTriangleValues(ATA),linspace(-500,500,100));
        xlim([-500 500])
        ylim([0 12000])
        title('Filter cov');        
        subplot(4,nStep,nStep + i);
        hist(upperTriangleValues(iATA),linspace(-0.01,0.01,100));  
        xlim([-0.01 0.01])
        ylim([0 5000])
        title('inv Filter cov');
        subplot(4,nStep,2*nStep + i);
        hist(diag(ATA),linspace(0,1500,100));
        xlim([0 1500])
        ylim([0 120])
        title('Filter var');
        subplot(4,nStep,3*nStep + i);
        hist(diag(iATA),linspace(0,0.03,100));
        xlim([0 0.03])
        ylim([0 120])
        title('inv Filter var');
        
        for j = 1:nFilt                
            OSIs(i,j) = orientationSelectivityIndex(tuning_curves(i,j,:),FRnonlin);
        end
%         hist(prefThetas(i,:),linspace(0,180,100));
%         title(sprintf('EM step %d',em_steps(i)))
%         pause
    end
    
    figure;
    errorbar(em_steps,mean(OSIs,2),std(OSIs,0,2)/sqrt(nFilt),'LineWidth',2);
    [s,p] = ttest(OSIs(1,:),OSIs(end,:))
    xlim([0 em_steps(nStep)+1]);
    xlabel('EM step #','FontSize',16);
    ylabel('OSI mean and s.e.m.','FontSize',16);
    figure;
    errorbar(em_steps,mean(postVars,2),std(postVars,0,2)/sqrt(nFilt),'LineWidth',2);
    xlim([0 em_steps(nStep)+1]);
    xlabel('EM step #','FontSize',16);
    ylabel('Response variance, mean and s.e.m.','FontSize',16);
    figure;
    errorbar(em_steps,mean(zmeans,2),std(zmeans,0,2)/sqrt(nFilt),'LineWidth',2);
    ylabel('Z MEAN mean and s.e.m.','FontSize',16);
     figure;
    errorbar(em_steps,mean(zvars,2),std(zvars,0,2)/sqrt(nFilt),'LineWidth',2);
    ylabel('Z VAR mean and s.e.m.','FontSize',16);
end
        