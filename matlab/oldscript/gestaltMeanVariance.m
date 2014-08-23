function gestaltMeanVariance()
    Dx = 64;
    k = 2;
    delta_x = 1;
    cc = gestaltCovariances(k,Dx,Dx);
    A = eye(Dx);
    ps = zeros(sqrt(Dx));
    ps(2:4,3) = 1;
    x = ps(:);
    g1 = 0:0.1:1;
    pixel_index = 2 * sqrt(Dx) + 6;
    other_index = 5 * sqrt(Dx) + 6; % this is inside another gestalt
    
    means = zeros(size(g1));
    variances = zeros(size(g1));
    
    other_means = zeros(size(g1));
    other_variances = zeros(size(g1));
    
    for i=1:size(g1,2)
        g = [g1(i);1-g1(i)];
        Cv = componentSum(g,cc);        
        pv_cov = inv((1/delta_x) * (A'*A) + inv(Cv));
        pv_mean = (1/delta_x) * pv_cov * A' * x;
        means(1,i) = pv_mean(pixel_index);
        variances(1,i) = pv_cov(pixel_index,pixel_index);
        other_means(1,i) = pv_mean(other_index);
        other_variances(1,i) = pv_cov(other_index,other_index);
    end
    
    subplot(1,2,1);
    plot(g1,means);
    hold on;
    plot(g1,other_means,'r');
    legend({'gestalt','other'});
    xlabel('Coefficient of gestalt covariance component');
    title('Means of V1 variables');
    
    subplot(1,2,2);
    plot(g1,variances);
    hold on;
    plot(g1,other_variances,'r');
    legend({'gestalt','other'});
    xlabel('Coefficient of gestalt covariance component');
    title('Variances of V1 variables');
end
    