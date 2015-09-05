function gsmOrientSelEM(iterfile,loadStuff,posterior,thetaRes)
    load(iterfile);
    nStep = length(A_iter);
    nFilt = size(A_iter{1},2);
    if loadStuff
        load('save_emorient.mat');
    else
        OSIs = zeros(nStep,nFilt);
        for i=1:nStep        
            [~,tuningCurves] = gsmOrientationSelectivity(A_iter{i},C,sigma_iter(i),thetaRes,false,'grating',~posterior,1,false);
            for j = 1:nFilt
                OSIs(i,j) = orientationSelectivityIndex(tuningCurves(j,:));
            end
        end
        save('bin/save_emorient.mat','OSIs');
    end
    errorbar(em_steps,mean(OSIs,2),std(OSIs,0,2),'LineWidth',2);
    xlim([0 em_steps(nStep)+20]);
    xlabel('EM step #','FontSize',16);
    ylabel('OSI mean and st.d.','FontSize',16);
end
        