function [means,stds,legends] = plotLearningCurves(names)
    means = [];
    stds = [];
    legends = {};
    for n=1:size(names,2)
        existing_param = true;
        act_param = 1;
        while existing_param
            runfiles = dir(sprintf('%s_iter_param%d_run*.mat',names{n},act_param));
            if isempty(runfiles)
                existing_param = false;
                continue;
            end            
            differences = [];
            for f = 1:size(runfiles,1)
                load(sprintf('%s_iter_param%d_run%d.mat',names{n},act_param,f));
                stepnum = size(state_sequence,2);
                actdiff = zeros(1,stepnum);
                for s = 1:stepnum
                    actdiff(1,s) = state_sequence{s}.difference_to_truth;
                end
                differences = [differences; actdiff];
            end
            means = [means; mean(differences)]; % this will break if step numbers differ
            stds = [stds; std(differences)]; % this will break if step numbers differ
            legends{size(legends,2)+1} = sprintf('%s_param%d',names{n},act_param);
            
            act_param = act_param + 1;
        end
    end
    X = repmat(1:size(means,2),size(means,1),1);
    errorbar(X,means,stds);
    %errorbar(X',log(means'),log(stds'));
    legend(legends);
end