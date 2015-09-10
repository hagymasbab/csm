%cd /home/umat/ogergo/a_learn_data
cd /home/banmi/gsm/run2
files = dir('*.mat');
A_iter = {};
sigma_iter = [];
em_steps = [];
for i=1:length(files)
    actFile = files(i);
    load(actFile.name);
    splitName = strsplit(actFile.name,'_');
    stepNum = str2num(splitName{end}(3:end-4));
%     if strcmp(splitName{end-2},'restart2k4')
%         stepNum = stepNum + 2400;
%     end
    if i ==1
        C = rho;
    end
    
    lastonly = false;
    done = false;
    if lastonly
        j = length(A_seq);        
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
    else
        j = 1;
        while ~done
            if isempty(A_seq{j})
                done = true;
            else
                A_iter{end+1} = A_seq{j};
                sigma_iter = [sigma_iter sigmaX_seq(j)];
                if isempty(em_steps)
                    em_steps = [em_steps 1];
                else
                    em_steps = [em_steps em_steps(end)+1];
                end
                j = j+1;
            end
        end
    end
end

[em_steps,sortIdx] = sort(em_steps);
sigma_iter = sigma_iter(sortIdx);
A_iter = A_iter(sortIdx);
        
cd /home/banmi/csm/matlab
save('bin/A_gsm_iter3.mat','C','A_iter','sigma_iter','em_steps');
