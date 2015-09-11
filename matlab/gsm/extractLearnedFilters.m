%cd /home/umat/ogergo/a_learn_data
cd /home/banmi/gsm/run2
files = dir('*.mat');

cycles_steps = [1 3 4 5 6; 10 8 20 25 30];

A_iter = {};
sigma_iter = [];
em_steps = [];
nFiles = size(cycles_steps,2);
for i=1:nFiles
    actname = sprintf('natim_learnA_2_dx256_du248_dz50_normem0_restart2k4_cycle%d_em%d.mat',cycles_steps(1,i),cycles_steps(2,i))
    %actFile = files(i);
    %load(actFile.name);
    load(actname)
    %splitName = strsplit(actFile.name,'_');
    splitName = strsplit(actname,'_');
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
save('bin/A_gsm_iter4.mat','C','A_iter','sigma_iter','em_steps');
