function plotLearning(benchmarkName,goal,matchComponents)
          
    summedGoal = componentSum(1,goal);
    Dv = size(goal{1},1);
    k = size(goal,2);
    
    if matchComponents
        % TODO
    end    
    
    firstruns = dir(sprintf('%s_iter_param*_run1.mat',benchmarkName));
    parNum = size(firstruns,1);
    runfiles = dir(sprintf('%s_iter_param1_run*.mat',benchmarkName));  
    runNum = size(runfiles,1);
    
    summedGoalParams = repmat(summedGoal(:),runNum,1);
    
    colNum = 1; % how many things will we plot for one parametrisation
            
    for act_param = 1:parNum
        sumInitParams = zeros(Dv*Dv*runNum,1);
        sumLearnedParams = zeros(Dv*Dv*runNum,1);
        for f = 1:runNum
            filename = sprintf('%s_iter_param%d_run%d.mat',benchmarkName,act_param,f);
            load(filename);
            initialComps = state_sequence{1}.estimated_components;
            learnedComps = state_sequence{end}.estimated_components;
            summedInitial = componentSum(1,initialComps);
            summedLearned = componentSum(1,learnedComps);
            sumInitParams((f-1)*Dv^2+1:f*Dv^2) = summedInitial(:);
            sumLearnedParams((f-1)*Dv^2+1:f*Dv^2) = summedLearned(:);
            if matchComponents
                % TODO
            end                        
        end
        
        subplot(parNum,colNum,(act_param-1)*colNum+1);
        scatter(summedGoalParams,sumLearnedParams,'bo');
        hold on;
        scatter(summedGoalParams,sumInitParams,'gx');        
        plot(xlim,xlim,'r-');
        if act_param == 1
            title('Summed parameters','FontSize',20);
            legend({'learned','initial'},'FontSize',16);
        end
        xlabel('True value','FontSize',16);
        ylabel('Estimated value','FontSize',16);
        
    end
end
%     
%     for i=1:5
%         filename = sprintf('test_iter_param1_run%d.mat',i);
%         load(filename);
%         cc1 = state_sequence{1}.estimated_components{1};
%         cc2 = state_sequence{1}.estimated_components{2};
%         if i==1 || i==3 || i==4
%             temp = cc1;
%             cc1 = cc2;
%             cc2 = temp;
%         end
%         sumc = componentSum(ones(ge1.k,1),{cc1,cc2});
%         starts_all1 = [starts_all1 cc1(:)' cc2(:)'];
%         starts_sum1 = [starts_sum1 sumc(:)'];
%         
%         cc1 = state_sequence{8}.estimated_components{1};
%         cc2 = state_sequence{8}.estimated_components{2};
%         if i==1 || i==3 || i==4
%             temp = cc1;
%             cc1 = cc2;
%             cc2 = temp;
%         end
%         sumc = componentSum(ones(ge1.k,1),{cc1,cc2});
%         learned_all1 = [learned_all1 cc1(:)' cc2(:)'];
%         learned_sum1 = [learned_sum1 sumc(:)'];
%                 
%         % -----------------------------------------------------------
%         
%         filename = sprintf('pret_iter_param1_run%d.mat',i);
%         load(filename);
%         cc1 = state_sequence{1}.estimated_components{1};
%         cc2 = state_sequence{1}.estimated_components{2};
%         if i==1 || i==3 || i==4
%             temp = cc1;
%             cc1 = cc2;
%             cc2 = temp;
%         end
%         sumc = componentSum(ones(ge2.k,1),{cc1,cc2});
%         starts_all2 = [starts_all2 cc1(:)' cc2(:)'];
%         starts_sum2 = [starts_sum2 sumc(:)'];
%         
%         cc1 = state_sequence{6}.estimated_components{1};
%         cc2 = state_sequence{6}.estimated_components{2};
%         if i==1 || i==2
%             temp = cc1;
%             cc1 = cc2;
%             cc2 = temp;
%         end
%         sumc = componentSum(ones(ge2.k,1),{cc1,cc2});
%         learned_all2 = [learned_all2 cc1(:)' cc2(:)'];
%         learned_sum2 = [learned_sum2 sumc(:)'];
%     end
%        
%     subplot(2,3,1);
%     scatter(goals_all1,learned_all1);    
%     hold on;
%     scatter(goals_all1,starts_all1,'x');    
%     
%     subplot(2,3,2);
%     scatter(goals_sum1,learned_sum1);
%     hold on;
%     scatter(goals_sum1,starts_sum1,'x');
%     
%     subplot(2,3,4);
%     scatter(goals_all2,learned_all2);    
%     hold on;
%     scatter(goals_all2,starts_all2,'x');    
%     
%     subplot(2,3,5);
%     scatter(goals_sum2,learned_sum2);
%     hold on;
%     scatter(goals_sum2,starts_sum2,'x');
%     
%     subplot(2,3,3);
%     plotDifferences('test',5,1);
%     yl = ylim;
%     subplot(2,3,6);
%     plotDifferences('pret',5,1);
%     ylim(yl);
% end