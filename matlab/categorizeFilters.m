function cat = categorizeFilters(stimulus,filters,components,model)
    cat.single_criteria = {'ovl','stim'};
    if strcmp(model,'csm')
        cat.pairwise_criteria = {'var','covar'};
    elseif strcmp(model,'csm')
        cat.pairwise_criteria = {'mean'};
    end
    
    nCrit = length(cat.criteria);
    
    % produce positive and negative criteria-fulfilling sets
    cat.criteria_positive = [];
    cat.criteria_negative = [];
    
    for cr = 1:nCrit
    end
    
    % produce all possible pairwise criteria from single ones
    
    % produce all combinations of binary criteria to get categories
    
    
    
    % calculate indices for category pairs
end

function res = ovl()