function gestaltPlotZCond(ge,n,V)
    zxmax = 7;
    zx = 0:0.1:zxmax;
    zlp = zeros(size(zx));
    for i = 1:size(zx,2)
        zlp(1,i) = gestaltLogPostZ(zx(1,i),n,V,ge);
    end
    plot(zx,exp(zlp));
    xlim([0 zxmax]);
    hold on;    
    plot([ge.Z(n,1); ge.Z(n,1)],ylim(),'r-');
    title(sprintf('Z = %f',ge.Z(n,1)))
end