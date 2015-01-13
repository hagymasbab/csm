function [mindiff,minperm,maxel_diff] = covcompRootMeanSquare(cc1,cc2,minperm,varargin)
    parser = inputParser;
    addParamValue(parser,'verbose',false,@islogical);
    parse(parser,varargin{:});        
    params = parser.Results;     
    % we should take the minimum over all possible permutations    
    if ~iscell(cc1)
        cc1 = {cc1};
        cc2 = {cc2};
    end
    k = size(cc1,2);
    d = size(cc1{1},1);
    if isempty(minperm);
        permutations = perms(1:k);
    else
        permutations = minperm;
    end
    mindiff = Inf;
    maxel_diff = Inf;
    for p=1:size(permutations,1)
        if params.verbose
            printCounter(p,'stringVal','RMS permuataion','maxVal',size(permutations,1),'newLine',true);
        end
        si = permutations(p,:);
        diff = 0;
        act_maxel_diff = 0;
        for j=1:k            
            other_index = si(1,j);
%             size(cc2{j})
%             size(cc1{other_index})
            squared_diffs = (cc2{j}-cc1{other_index})^2;
            act_maxel_diff = max([act_maxel_diff max(squared_diffs(:))]);
            
            diff = diff + sum(sum(squared_diffs));           
        end
        diff = diff / (k*d^2);
        if diff < mindiff
            mindiff = diff;
            maxel_diff = act_maxel_diff;
            minperm = si;
        end
    end
end