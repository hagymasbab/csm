function temporalSequence(ge,halfstim)
    randstim = randn(ge.B,ge.Dx);
    
    ge.obsVar = 1;
    
    ps = zeros(sqrt(ge.Dx));
    ps(2:6,3) = 1;
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
    v_stim = mean(v_samp(stim_idx));
    v_gest = mean(v_samp(gest_idx));
    v_conf = mean(v_samp(conf_idx));
    v_other = mean(v_samp(other_idx));
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
    v_stim = mean(v_samp(stim_idx));
    v_gest = mean(v_samp(gest_idx));
    v_conf = mean(v_samp(conf_idx));
    v_other = mean(v_samp(other_idx));
    v_act = [v_stim;v_gest;v_conf;v_other];
    v_means = [v_means v_act];
    
    % T4
    % update g
    
    plot([v_means;g_seq(1,:)]');
    h=legend('stimulated V1','gestalt V1','confound V1','other V1','higher');
    set(h,'Location','NorthEastOutside');
end    