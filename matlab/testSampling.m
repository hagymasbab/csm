function all_G = testSampling(ge,L,sampler,prior,plotLevel)
    right = 0;
    all_G = zeros(ge.N,L,ge.k);
    rightindices = [];
    fprintf('Datum %d/',ge.N);
    for n=1:ge.N
        printCounter(n);
        if plotLevel > 1
            if n==1
                sliceHandle = figure();
            else
                figure(sliceHandle);
            end
        end
        [s,~,zs] = gestaltGibbs(ge,n,L,'gSampler',sampler,'priorG',prior,'sampleRetry',100,'plotZ',plotLevel>1,'contrast',ge.contrast);
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
            posterior = @(g1,g2,g3) exp( gestaltLogPostG([g1;g2;g3;ge.G(n,4:end)'],ge.V(n,:,:),ge,prior,false) );
            
            if n==1
                ich = NaN;
            else
                ich = contourHandle;
            end
            [~,contourHandle] = plotFunctionProjections(posterior,0.1,ich);
            
            figure(contourHandle);
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
                zHandle = figure();
            else
                figure(zHandle);
            end
            clf;

            subplot(1,3,3);
            Vs = squeeze(mean(reshape(s(:,ge.k+1:end),L,ge.B,ge.Dv),2));
            Vtrue = squeeze(mean(ge.V(n,:,:),2));
            diff_mean = zeros(L,1);
            deviation = zeros(L,1);
            for i=1:L
                deviation(i,1) = mean(var(Vs(1:i,:)));
                rolling_mean = mean(Vs(1:i,:));
                diff_mean(i,1) = sqrt(sum((rolling_mean - Vtrue').^2));
            end
            plot([diff_mean deviation]);
            legend({'RMS of mean','mean of variance'});
            title('V sample convergence');           
                
            if ge.contrast       
                subplot(1,3,1);
                gestaltPlotZCond(ge,n,ge.V(n,:,:));
                hold on;
                ymax = ylim();
                ymax = ymax(2);
                scatter(zs,zeros(L,1) + ymax/2);
                
                subplot(1,3,2);
                plot(zs);
                hold on;
                plot(xlim(),[ge.Z(n,1); ge.Z(n,1)],'r-');
                ylim([0,2*ge.Z(n,1)]);                                                     
            end            
            
            pause;
        end
        
    end
    fprintf('\n%d/%d\n',right,ge.N);
    rightindices
    
    
end