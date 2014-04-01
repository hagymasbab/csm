function lp = gestaltPlotCondPostG(ge,V)
    % k=2
    gx = 0.01:0.01:0.99;
    lp = zeros(1,size(gx,2));
    for g=1:size(gx,2)
        lp(1,g) = gestaltLogPostG(gx(g),V,ge);
    end
    clf;
    plot(gx,exp(lp));
    %plot(gx,lp);
end