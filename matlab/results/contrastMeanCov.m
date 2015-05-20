function contrastMeanCov()
    Dv = 256;
    imdim = sqrt(Dv);
    % load filter set
    load(sprintf('filters_gabor_4or_%d.mat',Dv));
    
    % define covariance components
    cc{1} = eye(Dv);
    %filter_ids = [2*imdim+2*4 2*imdim+3*4 2*imdim+4*4];
    filter_ids = [146-4*16+8 146 146+4*16-8];
    cc{2} = 0.01 * eye(Dv);
    for f = 1:length(filter_ids)
        for ff = f+1:length(filter_ids)
            i = filter_ids(f);
            j = filter_ids(ff);
            cc{2}(i,j) = 1;
            cc{2}(j,i) = 1;
            cc{2}(i,i) = cc{2}(i,i) + 1;
            cc{2}(j,j) = cc{2}(j,j) + 1;
        end        
    end    
    %viewImageSet(cc);
    
    % TODO define mean components
    B = randn(Dv,2);    
    
    % construct base stimulus
    x_base = zeros(Dv,1);
    for f = 1:length(filter_ids)
        x_base = x_base + A(:,filter_ids(f));
    end
    %viewImage(x_base);
    
    % create contrast-adjusted stimuli
    contrasts = [0.01 1 10];
    
    x_rms = zeros(length(contrasts),1);
    corr_c = zeros(length(contrasts),1);
    corr_m = zeros(length(contrasts),1);
    cell1 = filter_ids(1);
    cell2 = filter_ids(2);
    
    for c = 1:length(contrasts)
        x_act = contrasts(c) * x_base;
        x_rms(c) = std(x_act(:)); 
        % get posterior covariances for each stimuli    
        [covc,covm] = posteriorCovariances(x_act,A,0.5,0.5,2,2,2,2,cc,B,200,'last');
        % get correlations and variances
        corrc = corrcov(covc);
        corrm = corrcov(covm);
        corr_c(c) = corrc(cell1,cell2);
        corr_m(c) = corrm(cell1,cell2);
    end
    % plot results
    plot(x_rms,[corr_c corr_m]);
end