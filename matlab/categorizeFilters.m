function cat = categorizeFilters(stimulus,filters,components,model)
    Dv = size(A,2);
    cat.criteria = {'ovl','stim','priorcov'};    
    
    nCrit = length(cat.criteria);
    
    cat.criteria_assignments = cell{1,nCrit};
    for cr = 1:nCrit
        positive = [];
        negative = [];
        for i=1:Dv
            for j=i+1:Dv
                assignment = feval(cat.criteria{cr},stimulus,filters,components,model);
                if assignment == 1
                    positive = [positive; i j];
                elseif assignment == -1
                    negative = [negative; i j];
                end
            end
        end
        crit_assign.positive = positive;
        crit_assign.negative = negative;
        cat.criteria_asignments{cr} = crit_assign;
    end
    
    
    % TODO produce all combinations of binary criteria to get categories
    
        
    % TODO calculate indices for category pairs
end

function res = ovl(stimulus,filters,components,model)
    % overlap of the two filters    
    f1 = filters(:,i1);
    f2 = filters(:,i2);
    overlap = abs(f1' * f2) / length(f1);
    % TODO figure out thresholds
    if overlap > 0.5
        res = 1;
    elseif overlap < 0.1
        res = -1;
    else
        res = 0;
    end
end

function res = stim(stimulus,filters,components,model)
    % similarity in activation by the stimulus

end

function res = priorcov(stimulus,filters,components,model)
    % absolute value of correlation implied by the prior
end