function temporalActivations()
    ge = gestaltCreate('nyolc','Dx',64,'filters','eye','B',10,'obsVar',0.1,'contrast',false);    
    ge.obsVar = 1;
    halfstim = gestaltStimulus(ge.Dx,ge.B,true,true);
    ge.X(1,:,:) = rand(ge.B,ge.Dx);
    ge.X(2,:,:) = halfstim;    
    ge.X(3,:,:) = halfstim;
    ge.X(4,:,:) = halfstim;
    
    ps = zeros(sqrt(ge.Dx));
    ps(2:4,3) = 1;
    ps = ps(:)';
    
    fullStim = zeros(sqrt(ge.Dv));
    fullStim(2:7,3) = 1;
    fullStim = fullStim(:)';
    
    gestaltBound = logical(fullStim - ps);
    other = logical(ones(1,ge.Dv) - fullStim);
    part = logical(ps);
    
    sampsteps = 1;
    nSamp = 4;
    nRestarts = 10;
    avg_stim = zeros(1,nSamp*2);
    avg_gest = zeros(1,nSamp*2);
    avg_other = zeros(1,nSamp*2);
    avg_g = zeros(1,nSamp*2);
    
    act_g = (1/ge.k) * ones(ge.k,1);
    %act_g = abs(rand(ge.k,1));
    for rest=1:nRestarts        
        for samp=1:nSamp          
            %fprintf('Samp %d\n',samp);
            s = gestaltGibbs(ge,samp,sampsteps,'initG',act_g,'priorG','gamma','verbose',0);
            act_g = s(sampsteps,1:ge.k)';
            % TEST
%             if samp < 3
%                 act_g = (1/ge.k) * ones(ge.k,1);
%             else
%                 act_g = [1; zeros(ge.k-1,1)];
%             end
            g = act_g(1);
%             while samp~=2 || g < 0.5
%                 s = gestaltGibbs(ge,samp,1);
%             end
%             g = s(1,1);
            v = reshape(s(sampsteps,ge.k+1:ge.k+ ge.B * ge.Dv),ge.B,ge.Dv);
            v_mean = mean(v);
            idx = zeros(1,nSamp*2);
            idx(1,2 * samp - 1:2 * samp) = 1;
            idx = logical(idx);            
            avg_stim(idx) = avg_stim(idx) + mean(v_mean(part));            
            avg_gest(idx) = avg_gest(idx) + mean(v_mean(gestaltBound));
            avg_other(idx) = avg_other(idx) + mean(v_mean(other));
            avg_g(idx) = avg_g(idx) + g;
        end                    
    end
    
    avg_stim = avg_stim / nRestarts;
    avg_gest = avg_gest / nRestarts;
    avg_other = avg_other / nRestarts;
    %avg_g = avg_g / nRestarts;
    avg_g = [avg_g(1,2:nSamp*2) avg_g(1,nSamp*2)] / nRestarts;    
    
    plot([avg_stim;avg_gest;avg_other;avg_g]');
    legend('stimulated V1','gestalt V1','other V1','higher');
end