function all_G = testSampling(ge,L,gsampler,zsampler,vsampler,plotLevel)
    if ge.k >= 3 && plotLevel > 0
        close all;
    end
    right = 0;
    all_G = zeros(ge.N,L,ge.k);
    rightindices = [];
    for n=1:ge.N
        %printCounter(n,'stringVal','Datum','maxVal',ge.N,'newLine',true);
        fprintf('Datum %d/%d\n',ge.N,n);
        if plotLevel > 1
            if n==1
                sliceHandle = figure();
            else
                figure(sliceHandle);
            end
        end
        [s,~,zs] = gestaltGibbs(ge,n,L,'gSampler',gsampler,'zSampler',zsampler,'vSampler',vsampler,'priorG','gamma', ...
            'sampleRetry',100,'plotZ',plotLevel>1,'contrast',ge.contrast,'verbose',1,'stepsize',0.2);
        all_G(n,:,:) = s(:,1:ge.k);
        meanG = mean(s(:,1:ge.k));
        [~,sampind] = max(meanG);
        [~,realind] = max(ge.G(n,:));
        if sampind == realind
            right = right + 1;
            rightindices = [rightindices n];
        else
            %ge.G(n,:)
        end
        
        if ge.k >= 3 && plotLevel > 0
            
            % we condition the posterior over g to the truth for the sake of plotting            
            posterior = @(g1,g2,g3) exp( gestaltLogPostG([g1;g2;g3;ge.G(n,4:end)'],ge.V(n,:,:),ge,'gamma',false) );
            
            if n==1
                contourHandle = figure('Units','normalized','OuterPosition',[0.1 0.1 0.6 0.4]);
            else
                figure(contourHandle);
            end
            clf;
            [~,contourHandle] = plotFunctionProjections(posterior,0.1,contourHandle);
            
            %figure(contourHandle);
            subplot(1,3,1);
            hold on;
            scatter(s(:,1),s(:,2));
            xlabel('comp1');
            ylabel('comp2');
            scatter(ge.G(n,1),ge.G(n,2),140,'rx','LineWidth',3);
            
            subplot(1,3,2);
            hold on;
            scatter(s(:,1),s(:,3));
            xlabel('comp1');
            ylabel('comp3');
            scatter(ge.G(n,1),ge.G(n,3),140,'rx','LineWidth',3);
            
            subplot(1,3,3);
            hold on;
            scatter(s(:,2),s(:,3));
            xlabel('comp2');
            ylabel('comp3');
            scatter(ge.G(n,2),ge.G(n,3),140,'rx','LineWidth',3);
            
            
            if n==1
                zHandle = figure('Units','normalized','OuterPosition',[0.1 0.6 0.8 0.4]);
            else
                figure(zHandle);
            end
            clf;
            
            Vs = squeeze(mean(reshape(s(:,ge.k+1:end),L,ge.B,ge.Dv),2));
            Vtrue = squeeze(mean(ge.V(n,:,:),2));
            [~,Vpostmean] = gestaltPostVRnd(ge,n,ge.G(n,:)',ge.Z(n,1),false);
            postdiff = sqrt(sum((Vpostmean - Vtrue).^2));
            
            diff_mean = zeros(L,1);
            moving_mean = zeros(L,1);
            deviation = zeros(L,1);
            diff_G = zeros(L,1);
            moving_G = zeros(L,1);
            for i=1:L
                deviation(i,1) = mean(var(Vs(1:i,:)));
                rolling_mean = mean(Vs(1:i,:));
                rolling_moving = mean(Vs(max(1,i-20):i,:));
                diff_mean(i,1) = sqrt(sum((rolling_mean - Vtrue').^2));
                moving_mean(i,1) = sqrt(sum((rolling_moving - Vtrue').^2));
                rolling_G = squeeze(mean(all_G(n,1:i,:),2));
                rolling_moving_G = squeeze(mean(all_G(n,max(1,i-20):i,:),2));
                %rolling_G - ge.G(n,:)'
                diff_G(i,1) = sqrt(sum((rolling_G - ge.G(n,:)').^2)');
                moving_G(i,1) = sqrt(sum((rolling_moving_G - ge.G(n,:)').^2));
            end
            subplot(1,4,3);
            plot([diff_mean moving_mean deviation]);
            hold on;
            plot(xlim(),[postdiff; postdiff],'r--');
            ylim([0,max(moving_mean)+2]);
            legend({'RMS of mean','RMS of rolling average','mean of variance'});
            title('V sample convergence');  
            subplot(1,4,4);
            plot([diff_G moving_G]);
            ylim([0,ge.k]);
            legend({'RMS of mean','RMS of rolling average'});
            title('G sample convergence');  
                
            if ge.contrast       
                subplot(1,4,1);
                gestaltPlotZCond(ge,n,ge.V(n,:,:));
                hold on;
                ymax = ylim();
                ymax = ymax(2);
                scatter(zs,zeros(L,1) + ymax/2);
                
                subplot(1,4,2);
                plot(zs);
                hold on;
                plot(xlim(),[ge.Z(n,1); ge.Z(n,1)],'r-');
                ylim([0,max(ge.Z(:,1))+0.1]);                                                     
            end            
            
            pause;
        end
        
    end
    fprintf('\n%d/%d\n',right,ge.N);
    rightindices
    
    
end