cd /home/umat/ogergo/a_learn_data
files = dir('*.mat');
A_iter = {};
sigma_iter = [];
em_steps = [];
for i=1:length(files)
    actFile = files(i);
    load(actFile.name);
    splitName = strsplit(actFile.name,'_');
    stepNum = str2num(splitName{end}(3:end-4));
    if strcmp(splitName{end-2},'restart2k4')
        stepNum = stepNum + 2400;
    end
    if i ==1
        C = rho;
    end
    j = length(A_seq);
    done = false;
    while ~done
        if isempty(A_seq{j})
            j = j-1;
        else
            done = true;
        end
    end
    A_iter{end+1} = A_seq{j};
    sigma_iter = [sigma_iter sigmaX_seq(j)];
    em_steps = [em_steps stepNum];
end

[em_steps,sortIdx] = sort(em_steps);
sigma_iter = sigma_iter(sortIdx);
A_iter = A_iter(sortIdx);
        
cd /home/banmi/csm/matlab
save('bin/A_gsm_iter2.mat','C','A_iter','sigma_iter','em_steps');
