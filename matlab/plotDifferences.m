function plotDifferences(name,nRun,nParam)
    for param = 1:nParam
        diff = [];
        norms = [];
        for n = 1:nRun
            actdiff = [];
            actnorms = [];
            filename = sprintf('%s_iter_param%d_run%d.mat',name,param,n);
            load(filename);    
            
            statenum = size(state_sequence,2);
            
            for state=1:statenum
                actdiff = [actdiff state_sequence{state}.difference_to_truth];
                actnorms = [actnorms cell2mat(state_sequence{state}.matrix_norms)];
            end

            diff = [diff;actdiff];
            norms = [norms;actnorms];

        end
        
        errorbar(mean(diff),std(diff))
        hold on;
    end
end