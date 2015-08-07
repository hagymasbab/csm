function cat = categorizeFilters(stimulus,filters,components,model,percentage,plotPairs)
    Dv = size(filters,2);
    cat.criteria = {'ovl','stim','prior'};    
    %cat.criteria = {'ovl','stim'};    
    %cat.criteria = {'prior'};    
    
    nCrit = length(cat.criteria);
    nExts = floor((Dv*(Dv-1)/2) * percentage);
    
    cat.criteria_assignments = cell(1,nCrit);
    for cr = 1:nCrit
        values = [];
        for i=1:Dv
            for j=i+1:Dv
                value = feval(cat.criteria{cr},stimulus,filters,components,model,i,j);
                values = [values; value i j];
            end
        end
        
        [~,topIndices,~] = maxNElements(values(:,1),nExts,[]);
        [~,bottomIndices,~] = maxNElements(-values(:,1),nExts,[]);
        positive = values(topIndices,2:3);
        negative = values(bottomIndices,2:3);
        
        crit_assign.positive = positive;
        crit_assign.negative = negative;
        cat.criteria_assignments{cr} = crit_assign;
    end
        
    % produce all combinations of binary criteria to get pair categories
    cat.categories = {};
    cat.category_assignments = {};
    % combinations = dec2bin(0:2^(length(cat.criteria))-1);
    combs = combinations([0 -1 1],length(cat.criteria));
    combs = combs(2:end,:); % drop the all zero line
    for c=1:size(combs,1)
        catname = [];
        filter_pairs = [];
        firstLoad = true;
        for cr=1:size(combs,2);            
            %if strcmp(combinations(c,cr),'0')
            if combs(c,cr) == -1
                catname = [catname cat.criteria{cr}(1)];            
                catname = [catname '-'];
                if firstLoad
                    filter_pairs = cat.criteria_assignments{cr}.negative;
                    firstLoad = false;
                else
                    filter_pairs = pair_intersection(filter_pairs,cat.criteria_assignments{cr}.negative);
                end
            elseif combs(c,cr) == 1
                catname = [catname cat.criteria{cr}(1)];            
                catname = [catname '+'];
                if firstLoad
                    filter_pairs = cat.criteria_assignments{cr}.positive;
                    firstLoad = false;
                else
                    filter_pairs = pair_intersection(filter_pairs,cat.criteria_assignments{cr}.positive);
                end
            end            
        end
        cat.categories{end+1} = catname;
        cat.category_assignments{end+1} = filter_pairs;
        if plotPairs
            nRow = size(combs,1);
            nCol = nRow;
            actCol = min(nCol,size(cat.category_assignments{c},1));
            for i = 1:actCol
                act_idx1 = cat.category_assignments{c}(i,1);
                act_idx2 = cat.category_assignments{c}(i,2);
                act_f1 = filters(:,act_idx1);
                act_f2 = filters(:,act_idx2);
                act_img = act_f1 + act_f2;
                subplot(nRow,nCol,(c-1)*nCol+i);
                viewImage(act_img,'useMax',true);      
                xlabel(sprintf('%d %d',act_idx1,act_idx2));
            end
            subplot(nRow,nCol,(c-1)*nCol+1);
            %viewImage(zeros(size(filters,1)));
            ylabel(cat.categories{c});
        end
    end
end
            

function overlap = ovl(stimulus,filters,components,model,i1,i2)
    % overlap of the two filters    
    f1 = filters(:,i1);
    f2 = filters(:,i2);
    overlap = abs(f1' * f2) / length(f1);
end

function diff = stim(stimulus,filters,components,model,i1,i2)
    % similarity in activation by the stimulus
    f1 = filters(:,i1);
    f2 = filters(:,i2);
    act1 = f1' * stimulus;
    act2 = f2' * stimulus;
    diff = abs(act1-act2);
end

function corr = prior(stimulus,filters,components,model,i1,i2)
    % absolute value of (some sort of) correlation implied by the prior
    if strcmp(model,'GSM')
        corr = components(i1,i2);
    elseif strcmp(model,'CSM')
        C = componentSum(1,components);
        corr = C(i1,i2);
    elseif strcmp(model,'MSM')
        B_f1 = components(i1,:);
        B_f2 = components(i2,:);
        corr = max(B_f1 .* B_f2); 
    else
        error('unsupported model type: %s',model);
    end        
end