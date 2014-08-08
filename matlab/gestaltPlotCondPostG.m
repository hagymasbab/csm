function lp = gestaltPlotCondPostG(ge,V,prior,precision)
    clf;
    % k=2
    gx = 0.01:0.01:0.99;
    gy = 0.01:0.01:0.99;
    lp = zeros(1,size(gx,2));
    lp2dim = zeros(size(gx,2),size(gy,2));
    for g=1:size(gx,2)
        if strcmp(prior,'dirichlet')
            lp(1,g) = gestaltLogPostG(gx(g),V,ge,prior,precision);
        else
            for gg=1:size(gy,2)
                lp2dim(g,gg) = gestaltLogPostG([gx(g);gy(gg)],V,ge,prior,precision);
            end
        end                
    end
    
    if ~strcmp(prior,'dirichlet')
        lp1 = sum(lp2dim,2);
        lp2 = sum(lp2dim,1);
        
    else
        lp1 = lp;
        lp2 = ones(1,size(gx,2)) - lp1; 
    end
    
    plot(gx,lp1);
    hold on;
    plot(gy,lp2,'r');
end