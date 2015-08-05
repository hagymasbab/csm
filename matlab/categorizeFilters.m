function cat = categorizeFilters(stimulus,filters,components,model)
    Dv = size(filters,2);
    %cat.criteria = {'ovl','stim','priorcov'};    
    cat.criteria = {'ovl','stim'};    
    
    nCrit = length(cat.criteria);
    
    cat.criteria_assignments = cell(1,nCrit);
    for cr = 1:nCrit
        positive = [];
        negative = [];
        res = [];
        for i=1:Dv
            for j=i+1:Dv
                [assignment,value] = feval(cat.criteria{cr},stimulus,filters,components,model,i,j);
                res = [res value];
                if assignment == 1
                    positive = [positive; i j];
                elseif assignment == -1
                    negative = [negative; i j];
                end
            end
        end
        crit_assign.positive = positive;
        crit_assign.negative = negative;
        cat.criteria_assignments{cr} = crit_assign;
        %hist(res,100);pause
    end
        
    % produce all combinations of binary criteria to get pair categories
    cat.categories = {};
    cat.category_assignments = {};
    combinations = dec2bin(0:2^(length(cat.criteria))-1);
    for c=1:size(combinations,1)
        catname = [];
        filter_pairs = [];
        for cr=1:size(combinations,2);
            catname = [catname cat.criteria{cr}];            
            if strcmp(combinations(c,cr),'0')
                catname = [catname '-'];
                if cr == 1
                    filter_pairs = cat.criteria_assignments{cr}.negative;
                else
                    filter_pairs = pair_intersection(filter_pairs,cat.criteria_assignments{cr}.negative);
                end
            else
                catname = [catname '+'];
                if cr == 1
                    filter_pairs = cat.criteria_assignments{cr}.positive;
                else
                    filter_pairs = pair_intersection(filter_pairs,cat.criteria_assignments{cr}.positive);
                end
            end            
        end
        cat.categories{end+1} = catname;
        cat.category_assignments{end+1} = filter_pairs;
    end
end


        
    

function [res,overlap] = ovl(stimulus,filters,components,model,i1,i2)
    % overlap of the two filters    
    f1 = filters(:,i1);
    f2 = filters(:,i2);
    overlap = abs(f1' * f2) / length(f1);
    % TODO figure out thresholds
    if overlap > 0.5e-3
        res = 1;
    elseif overlap < 0.05e-3
        res = -1;
    else
        res = 0;
    end
end

function [res,diff] = stim(stimulus,filters,components,model,i1,i2)
    % similarity in activation by the stimulus
    f1 = filters(:,i1);
    f2 = filters(:,i2);
    act1 = f1' * stimulus;
    act2 = f2' * stimulus;
    diff = abs(act1-act2);
    % TODO figure out thresholds
    if diff > 1
        res = 1;
    elseif diff < 0.1
        res = -1;
    else
        res = 0;
    end
end

function res = priorcov(stimulus,filters,components,model,i1,i2)
    % absolute value of correlation implied by the prior
end