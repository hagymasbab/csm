function all_G = testSampling(ge,L,sampler,prior,plotsamples)
    right = 0;
    all_G = zeros(ge.N,L,ge.k);
    rightindices = [];
    fprintf('Datum %d/',ge.N);
    for n=1:ge.N
        printCounter(n);
        [s,~,zs] = gestaltGibbs(ge,n,L,'gSampler',sampler,'priorG',prior,'sampleRetry',100,'plot',0,'contrast',ge.contrast);
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
        
        if ge.k >= 3 && plotsamples
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
            
            if ge.contrast
                if n==1
                    zHandle = figure();
                else
                    figure(zHandle);
                end
                
                clf;
                zxmax = 5;
                zx = 0:0.1:zxmax;
                zlp = zeros(size(zx));
                for i = 1:size(zx,2)
                    zlp(1,i) = gestaltLogPostZ(zx(1,i),n,ge.V(n,:,:),ge);
                end
                plot(zx,zlp);
                xlim([0 zxmax]);
                ymax = ylim();
                ymax = ymax(2);
                hold on;
                scatter(zs,zeros(L,1) + ymax/2);
                limits = ylim();
                plot([ge.Z(n,1); ge.Z(n,1)],limits,'r-');
                title(sprintf('Z = %f',ge.Z(n,1)))
            end            
            
            pause;
        end
        
    end
    fprintf('\n%d/%d\n',right,ge.N);
    rightindices
    
    
end