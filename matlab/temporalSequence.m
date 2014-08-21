function temporalSequence(ge,halfstim)
    randstim = randn(ge.B,ge.Dx);
    
    ge.obsVar = 1;
    sAA = (1/ge.obsVar) * ge.AA;
    
    ps = zeros(sqrt(ge.Dx));
    ps(2:5,3) = 1;
    ps = ps(:)';
    
    fullStim = zeros(sqrt(ge.Dv));
    fullStim(2:7,3) = 1;
    fullStim = fullStim(:)';
    
    confound = zeros(sqrt(ge.Dv));
    % TODO set confound
    confound = confound(:)';
    
    gest_idx = logical(fullStim - ps);
    stim_idx = logical(ps);
    conf_idx = logical(confound);
    other_idx = logical(ones(1,ge.Dv) - fullStim - confound);
    idxs = {stim_idx,gest_idx,conf_idx,other_idx};
    groupnum = 4;
    
    covidxs = cell{1,groupnum};
    for i=1:groupnum
        act = zeros(ge.Dv);
        for j=1:ge.Dv
            for k=j+1:ge.Dv
                if idxs{i}(j) == 1 && idxs{i}(k) ==1
                    act(j,k) = 1;
                end
            end
        end
        covidxs{i} = act;
    end
    
    
    v_means = []; % sim will be 4xT (stimulus,gestalt,confound,other)
    v_stds = [];
    g_seq = [];
    
    % T0
    % calculate baseline mean and variance of g from the prior
    % calculate mean and variance of v if there is no signal in the input
    g_init = (1/ge.k) * ones(ge.k,1);
    % TEST
    g_init = [0.1;0.9];
    g_seq = [g_seq g_init];
    v_means = [v_means zeros(4,1)];
    
    % T1
    % set the input to a partial activation of a gestalt plus unexplained
    % confounding variance
    % calculate mean and variance of v conditioned on the input and g
    g_seq = [g_seq g_init];
    ge.X(1,:,:) = halfstim;
    v_samp_batch = gestaltPostVRnd(ge,1,g_init,1,false);
    v_samp = mean(v_samp_batch);
%     v_stim = mean(v_samp(stim_idx));
%     v_gest = mean(v_samp(gest_idx));
%     v_conf = mean(v_samp(conf_idx));
%     v_other = mean(v_samp(other_idx));
    
    cv = componentSum(g_init,ge.cc);
    actcorr = corrcov(inv(sAA + inv(cv)));
    v_stim = meanCorr(v_samp_batch(:,stim_idx));
    v_gest = meanCorr(v_samp_batch(:,gest_idx));
    v_conf = meanCorr(v_samp_batch(:,conf_idx));
    v_other = meanCorr(v_samp_batch(:,other_idx));
    
    v_act = [v_stim;v_gest;v_conf;v_other];
    v_means = [v_means v_act];
    
    
    % T2
    % calculate mean and variance of g conditioned on v
    logpdf = @(g) gestaltLogPostG(g,v_samp_batch,ge,'dirichlet',false); 
    g_part = g_init(1:ge.k-1,1);
    for i=1:5
        [g_part,~] = sliceSample(g_part,logpdf,0.05,'limits',[0,1]);
    end
    g_act = [g_part; 1-sum(g_part)];
    % TEST
    %g_act = [1;0];
    g_seq = [g_seq g_act];
    v_means = [v_means v_act];
    
    % T3
    % set the input signal back to nothing
    % update v
    g_seq = [g_seq g_act];
    %ge.X(1,:,:) = randstim;
    ge.X(1,:,:) = halfstim;
    v_samp_batch = gestaltPostVRnd(ge,1,g_act,1,false);
    v_samp = mean(v_samp_batch);
%     v_stim = mean(v_samp(stim_idx));
%     v_gest = mean(v_samp(gest_idx));
%     v_conf = mean(v_samp(conf_idx));
%     v_other = mean(v_samp(other_idx));

    v_stim = meanCorr(v_samp_batch(:,stim_idx));
    v_gest = meanCorr(v_samp_batch(:,gest_idx));
    v_conf = meanCorr(v_samp_batch(:,conf_idx));
    v_other = meanCorr(v_samp_batch(:,other_idx));

    v_act = [v_stim;v_gest;v_conf;v_other];
    v_means = [v_means v_act];
    
    % T4
    % update g
    
    plot([v_means;g_seq(1,:)]');
    h=legend('stimulated V1','gestalt V1','confound V1','other V1','higher');
    set(h,'Location','NorthEastOutside');
end    

function mr = meanCorr(X)
%     co = corrcoef(X);
%     ut = triu(co) - diag(diag(co));
    mr = mean(X(:));
end