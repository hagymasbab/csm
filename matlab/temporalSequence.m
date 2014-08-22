function temporalSequence(ge,stim,quantity,remove)
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
    % set confound
    confound(2:5,4) = 1;
    confound = confound(:)';
    
    gest_idx = logical(fullStim - ps);
    stim_idx = logical(ps);
    conf_idx = logical(confound);
    other_idx = logical(ones(1,ge.Dv) - fullStim - confound);
    idxs = {stim_idx,gest_idx,conf_idx,other_idx};
    groupnum = 4;
    
    covidxs = cell(1,groupnum);
    for i=1:groupnum
        act = zeros(ge.Dv);
        for j=1:ge.Dv
            for k=j+1:ge.Dv
                if idxs{i}(j) == 1 && idxs{i}(k) ==1
                    act(j,k) = 1;
                end
            end
        end
        covidxs{i} = logical(act);
    end
    
    
    v_means = []; % sim will be 4xT (stimulus,gestalt,confound,other)
    v_stds = [];
    g_seq = [];
    
    % T0
    % calculate baseline mean and variance of g from the prior
    % calculate mean and variance of v if there is no signal in the input
    g_act = (1/ge.k) * ones(ge.k,1);
    % TEST
    g_act = [0.1;0.9];
    g_seq = [g_seq g_act];
    v_means = [v_means zeros(4,1)];
    
    % T1
    % set the input to a partial activation of a gestalt plus unexplained
    % confounding variance
    % calculate mean and variance of v conditioned on the input and g
    g_seq = [g_seq g_act];
    ge.X(1,:,:) = stim;
    v_samp_batch = gestaltPostVRnd(ge,1,g_act,1,false);
    v_samp = mean(v_samp_batch);
    
    if strcmp(quantity,'mean')
        v_stim = mean(v_samp(stim_idx));
        v_gest = mean(v_samp(gest_idx));
        v_conf = mean(v_samp(conf_idx));
        v_other = mean(v_samp(other_idx));
        v_act = [v_stim;v_gest;v_conf;v_other];
    elseif strcmp(quantity,'noisecorr')
        cv = componentSum(g_act,ge.cc);
        actcorr = corrcov(inv(sAA + inv(cv)));
        
        v_act = [];
        for i = 1:groupnum
            selected_corr = actcorr(covidxs{i});
            v_act = [v_act; mean(selected_corr(:))];
        end
    elseif strcmp(quantity,'allcorr')
        v_act = [];
        for i = 1:groupnum
            mc = meanCorr(v_samp_batch(:,idxs{i}));
            v_act = [v_act; mc];
        end
    end
    v_means = [v_means v_act];        
    
    % T2
    % calculate mean and variance of g conditioned on v
    logpdf = @(g) gestaltLogPostG(g,v_samp_batch,ge,'dirichlet',false); 
    g_part = g_act(1:ge.k-1,1);
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
    if remove
        ge.X(1,:,:) = randstim;
    else
        ge.X(1,:,:) = stim;
    end
    v_samp_batch = gestaltPostVRnd(ge,1,g_act,1,false);
    v_samp = mean(v_samp_batch);

    if strcmp(quantity,'mean')
        v_stim = mean(v_samp(stim_idx));
        v_gest = mean(v_samp(gest_idx));
        v_conf = mean(v_samp(conf_idx));
        v_other = mean(v_samp(other_idx));
        v_act = [v_stim;v_gest;v_conf;v_other];
    elseif strcmp(quantity,'noisecorr')
        cv = componentSum(g_act,ge.cc);
        actcorr = corrcov(inv(sAA + inv(cv)));
        
        v_act = [];
        for i = 1:groupnum
            selected_corr = actcorr(covidxs{i});
            v_act = [v_act; mean(selected_corr(:))];
        end
    elseif strcmp(quantity,'allcorr')
        v_act = [];
        for i = 1:groupnum
            mc = meanCorr(v_samp_batch(:,idxs{i}));
            v_act = [v_act; mc];
        end
    end
    v_means = [v_means v_act];   
    
    % T4
    % update g
    
    ph = plot([v_means;g_seq(1,:)]','LineWidth',1.5);
    set(ph(5),'linewidth',3);
    hold on;
    lh=legend('stimulated V1','gestalt V1','confound V1','other V1','higher');
    set(lh,'Location','NorthEastOutside');
    
    if remove
        areaxlim = [2 3];
    else
        areaxlim = [2 4];
    end
    yl = ylim;
    H = area([areaxlim fliplr(areaxlim)],[yl(2) yl(2) yl(1) yl(1)],'LineStyle','none');
    h=get(H,'children');
    set(h,'FaceAlpha',0.05,'FaceColor',[0 0 1]);
end


function mr = meanCorr(X)
    dim = size(X,1);
    co = corrcoef(X);
    ut = triu(co) - diag(diag(co));
    mr = mean(ut(:)) * (dim^2 / (dim^2 / 2 + dim / 2));
end