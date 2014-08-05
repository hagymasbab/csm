function temporalActivations()
    ge = gestaltCreate('nyolc','Dx',64,'filters','eye','B',10,'obsVar',0.1,'contrast',false);    
    ge.obsVar = 1;
    halfstim = gestaltStimulus(ge.Dx,ge.B,true,true);
    ge.X(1,:,:) = halfstim;
    
    ps = zeros(sqrt(ge.Dx));
    ps(2:4,3) = 1;
    ps = ps(:)';
    
    fullStim = zeros(sqrt(ge.Dv));
    fullStim(2:7,3) = 1;
    fullStim = fullStim(:)';
    
    gestaltBound = logical(fullStim - ps);
    other = logical(ones(1,ge.Dv) - fullStim);
    part = logical(ps);
    
    nSamp = 4;
    nRestarts = 10;
    avg_stim = zeros(1,nSamp);
    avg_gest = zeros(1,nSamp);
    avg_other = zeros(1,nSamp);
    avg_g = zeros(1,nSamp);
    for rest=1:nRestarts
        s = gestaltGibbs(ge,1,nSamp);
        for samp=1:nSamp
            v = reshape(s(samp,ge.k+1:ge.k+ ge.B * ge.Dv),ge.B,ge.Dv);
            v_mean = mean(v);
            avg_stim(1,samp) = avg_stim(1,samp) + mean(v_mean(part));
            avg_gest(1,samp) = avg_gest(1,samp) + mean(v_mean(gestaltBound));
            avg_other(1,samp) = avg_other(1,samp) + mean(v_mean(other));
            avg_g(1,samp) = avg_g(1,samp) + s(samp,1);
        end                    
    end
    
    avg_stim = avg_stim / nRestarts;
    avg_gest = avg_gest / nRestarts;
    avg_other = avg_other / nRestarts;
    avg_g = [0.5*nRestarts avg_g(1,2:nSamp)] / nRestarts;    
    
    plot([avg_stim;avg_gest;avg_other;avg_g]');
    legend('stimulated V1','gestalt V1','other V1','higher');
end