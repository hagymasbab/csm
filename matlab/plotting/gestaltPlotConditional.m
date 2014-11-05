function gestaltPlotConditional(g,idx,V,ge,prior,stepsize)
    gx = 0.01:stepsize:0.99;
    lp = zeros(1,size(gx,2));
    for gi=1:size(gx,2)
        lp(1,gi) = gestaltLogCondPostG(gx(gi),g,idx,V,ge,prior,false); 
    end
    
    plot(gx,exp(lp));
    xlim([-0.06 1]);
    title(sprintf('Component %d',idx));
end