function plotLearning(ge1,ge2)
    %pret: 1 2 - ford
    %test: 1 3 4 - ford        
    
    cv1 = componentSum(ones(ge1.k,1),ge1.cc);
    cv2 = componentSum(ones(ge2.k,1),ge2.cc);
    
    goals_all1 = repmat([ge1.cc{1}(:)' ge1.cc{1}(:)'],1,5);
    goals_all2 = repmat([ge2.cc{1}(:)' ge2.cc{1}(:)'],1,5);
    goals_sum1 = repmat(cv1(:)',1,5);
    goals_sum2 = repmat(cv2(:)',1,5);
    
    starts_all1 = [];
    starts_all2 = [];
    starts_sum1 = [];
    starts_sum2 = [];
    
    learned_all1 = [];
    learned_all2 = [];
    learned_sum1 = [];
    learned_sum2 = [];
    
    for i=1:5
        filename = sprintf('test_iter_param1_run%d.mat',i);
        load(filename);
        cc1 = state_sequence{1}.estimated_components{1};
        cc2 = state_sequence{1}.estimated_components{2};
        if i==1 || i==3 || i==4
            temp = cc1;
            cc1 = cc2;
            cc2 = temp;
        end
        sumc = componentSum(ones(ge1.k,1),{cc1,cc2});
        starts_all1 = [starts_all1 cc1(:)' cc2(:)'];
        starts_sum1 = [starts_sum1 sumc(:)'];
        
        cc1 = state_sequence{8}.estimated_components{1};
        cc2 = state_sequence{8}.estimated_components{2};
        if i==1 || i==3 || i==4
            temp = cc1;
            cc1 = cc2;
            cc2 = temp;
        end
        sumc = componentSum(ones(ge1.k,1),{cc1,cc2});
        learned_all1 = [learned_all1 cc1(:)' cc2(:)'];
        learned_sum1 = [learned_sum1 sumc(:)'];
                
        % -----------------------------------------------------------
        
        filename = sprintf('pret_iter_param1_run%d.mat',i);
        load(filename);
        cc1 = state_sequence{1}.estimated_components{1};
        cc2 = state_sequence{1}.estimated_components{2};
        if i==1 || i==3 || i==4
            temp = cc1;
            cc1 = cc2;
            cc2 = temp;
        end
        sumc = componentSum(ones(ge2.k,1),{cc1,cc2});
        starts_all2 = [starts_all2 cc1(:)' cc2(:)'];
        starts_sum2 = [starts_sum2 sumc(:)'];
        
        cc1 = state_sequence{6}.estimated_components{1};
        cc2 = state_sequence{6}.estimated_components{2};
        if i==1 || i==2
            temp = cc1;
            cc1 = cc2;
            cc2 = temp;
        end
        sumc = componentSum(ones(ge2.k,1),{cc1,cc2});
        learned_all2 = [learned_all2 cc1(:)' cc2(:)'];
        learned_sum2 = [learned_sum2 sumc(:)'];
    end
       
    subplot(2,3,1);
    scatter(goals_all1,learned_all1);    
    hold on;
    scatter(goals_all1,starts_all1,'x');    
    
    subplot(2,3,2);
    scatter(goals_sum1,learned_sum1);
    hold on;
    scatter(goals_sum1,starts_sum1,'x');
    
    subplot(2,3,4);
    scatter(goals_all2,learned_all2);    
    hold on;
    scatter(goals_all2,starts_all2,'x');    
    
    subplot(2,3,5);
    scatter(goals_sum2,learned_sum2);
    hold on;
    scatter(goals_sum2,starts_sum2,'x');
    
    subplot(2,3,3);
    plotDifferences('test',5,1);
    yl = ylim;
    subplot(2,3,6);
    plotDifferences('pret',5,1);
    ylim(yl);
end