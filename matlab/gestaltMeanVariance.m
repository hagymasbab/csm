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
    
    means = zeros(size(g1));
    variances = zeros(size(g1));
    
    for i=1:size(g1,2)
        g = [g1(i);1-g1(i)];
        Cv = componentSum(g,cc);        
        pv_cov = inv((1/delta_x) * (A'*A) + inv(Cv));
        pv_mean = (1/delta_x) * pv_cov * A' * x;
        means(1,i) = pv_mean(pixel_index);
        variances(1,i) = pv_cov(pixel_index,pixel_index);
    end
    
    plot(g1,means);
    hold on;
    plot(g1,variances,'r');
    legend({'mean','variance'});
    xlabel('Coefficient of relevant covariance component');
    title('Conditional density parameters of a gestalt-activated V1 variable');
end
    