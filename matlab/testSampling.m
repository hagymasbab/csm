function all_G = testSampling(ge,L,sampler,prior)
    right = 0;
    all_G = zeros(ge.N,L,ge.k);
    rightindices = [];
    fprintf('Datum %d/',ge.N);
    for n=1:ge.N
        printCounter(n);
        s = gestaltGibbs(ge,n,L,'gSampler',sampler,'priorG',prior,'sampleRetry',10,'plot',0);
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
        
        if ge.k == 3
            % we condition the posterior over g to the truth for the sake of plotting
            clf;
            posterior = @(g1,g2,g3) exp( gestaltLogPostG([g1;g2;g3],ge.V(n,:,:),ge,prior,false) );
            plotFunctionProjections(posterior,0.1);                        
            subplot(1,3,1);
            hold on;
            scatter(s(:,1),s(:,2));
            xlabel('comp1');
            ylabel('comp2');
            
            subplot(1,3,2);
            hold on;
            scatter(s(:,1),s(:,3));
            xlabel('comp1');
            ylabel('comp3');
            
            subplot(1,3,3);
            hold on;
            scatter(s(:,2),s(:,3));
            xlabel('comp2');
            ylabel('comp3');
            
            pause;
        end
        
    end
    fprintf('\n%d/%d\n',right,ge.N);
    rightindices
    
    
end