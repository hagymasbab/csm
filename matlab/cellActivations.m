function cellActivations(ge)
    ge.obsVar = 1;

    % create partial stimuli
    % there will be 20 of those, every possible way we can select 3
    % elements out of 6
    nStim = 1;    
    
    fullStim = zeros(sqrt(ge.Dv));
    if ge.k > 10
        fullStim(2:6,1) = 1;
    else
        fullStim(2:7,3) = 1;
    end
    fullStim = fullStim(:)';
    partStim = zeros(nStim,ge.Dx);        
    
    randPartStim = zeros(nStim,ge.B,ge.Dx);
    
    % draw samples for each, start over 10 times
    nSamp = 50;
    nRestarts = 100;
    samples = zeros(nStim,nRestarts,nSamp,ge.k + ge.B * ge.Dv);
    
    % average rates and variances of stimulated, gestalt-bound and
    % unrelated units, for all samples, and separately for correctly and
    % incorrectly leaning g-samples
    correct_num = 0;
    all_vmean = zeros(2,3);
    correct_vmean = zeros(2,3);
    incorr_vmean = zeros(2,3);
    all_vstds = zeros(2,3);
    correct_vstds = zeros(2,3);
    incorr_vstds = zeros(2,3);
    
    for stim=1:nStim
        
        % TODO do the rest of possible partial stimuli
        ps = zeros(sqrt(ge.Dx));
        if ge.k>10
            ps(3:5,1) = 1;
        else
            ps(2:4,3) = 1;
        end
        partStim(stim,:) = ps(:)';
        randPartStim(stim,:,:) = gestaltStimulus(ge.Dx,ge.B,true,true);
        ge.X(1,:,:) = randPartStim(stim,:,:);
        
        fprintf('Restart %d/',nRestarts);
        for samp=1:nRestarts
            printCounter(samp);
            initG = (1/ge.k) * ones(ge.k,1);
            valid = false;
            while ~valid
                [act_samp,rr] = gestaltGibbs(ge,1,nSamp,'initG',initG,'contrast',false);
                if rr ~= -1
                    valid = true;
                end
            end
            samples(stim,samp,:,:) = act_samp;
            
            g = reshape(act_samp(:,1:ge.k),nSamp,ge.k);
            v = reshape(act_samp(:,ge.k+1:ge.k+ ge.B * ge.Dv),nSamp*ge.B,ge.Dv);
            
            v_mean = mean(v);
            v_std = std(v);
            
            gestaltBound = logical(fullStim - partStim(stim));
            other = logical(ones(1,ge.Dv) - fullStim);
            part = logical(partStim(stim,:));
            
            stimulated_mean_mean = mean(v_mean(part));
            gestalt_mean_mean = mean(v_mean(gestaltBound));
            other_mean_mean = mean(v_mean(other));
            stimulated_mean_std = std(v_mean(part));
            gestalt_mean_std = std(v_mean(gestaltBound));
            other_mean_std = std(v_mean(other));
            vmeans = [stimulated_mean_mean gestalt_mean_mean other_mean_mean;stimulated_mean_std gestalt_mean_std other_mean_std];
            
            stimulated_std_mean = mean(v_std(part));
            gestalt_std_mean = mean(v_std(gestaltBound));
            other_std_mean = mean(v_std(other));
            stimulated_std_std = std(v_std(part));
            gestalt_std_std = std(v_std(gestaltBound));
            other_std_std = std(v_std(other));
            vstds = [stimulated_std_mean gestalt_std_mean other_std_mean;stimulated_std_std gestalt_std_std other_std_std];
            
            all_vmean = all_vmean + vmeans;
            all_vstds = all_vstds + vstds;
            
            [~,max_sample] = max(mean(g));
            if ge.k > 10
                max_truth = 3;
            else
                max_truth = 1;
            end
            if max_sample == max_truth
                correct_num = correct_num + 1;
                correct_vmean = correct_vmean + vmeans;
                correct_vstds = correct_vstds + vstds;
            else
                incorr_vmean = incorr_vmean + vmeans;
                incorr_vstds = incorr_vstds + vstds;
            end
        end
    end
    fprintf('\n');
    
    N = nStim*nRestarts;
    all_vmean = all_vmean / N;
    all_vstds = all_vstds / N;
    
    correct_vmean = correct_vmean / correct_num;
    correct_vstds = correct_vstds / correct_num;
        
    incorr_vmean = incorr_vmean / (N-correct_num);
    incorr_vstds = incorr_vstds / (N-correct_num);
    
    labels = {'stimulus','gestalt','other'};
    
    subplot(2,3,1);
    barwitherr(all_vmean(2,:),all_vmean(1,:));
    set(gca,'XTickLabel',labels);
    title(sprintf('All samples mean (%d)',N));
    
    subplot(2,3,4);
    barwitherr(all_vstds(2,:),all_vstds(1,:));
    set(gca,'XTickLabel',labels);
    title('All samples deviation');
    
    subplot(2,3,2);
    barwitherr(correct_vmean(2,:),correct_vmean(1,:));
    set(gca,'XTickLabel',labels);
    title(sprintf('Correct samples mean (%d)',correct_num));
    
    subplot(2,3,5);
    barwitherr(correct_vstds(2,:),correct_vstds(1,:));
    set(gca,'XTickLabel',labels);
    title('Correct samples deviation');
    
    subplot(2,3,3);
    barwitherr(incorr_vmean(2,:),incorr_vmean(1,:));
    set(gca,'XTickLabel',labels);
    title('Incorrect samples mean');
    
    subplot(2,3,6);
    barwitherr(incorr_vstds(2,:),incorr_vstds(1,:));
    set(gca,'XTickLabel',labels);
    title('Incorrect samples deviation');
    
end
    
    
    
    