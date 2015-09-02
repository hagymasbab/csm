function correlationPlot(x1,x2)
    cr = corrcoef(x1,x2);
    c = cr(1,2);
    def_colors = get(groot,'DefaultAxesColorOrder');
    reddish_color = def_colors(2,:);
    rfont = 20;
    p1 = polyfit(x1,x2,1);
    xlims = xlim();
    ylims = ylim();
    range = xlims(1):0.01:xlims(2);
    plot(range,polyval(p1,range),'LineWidth',3,'Color',reddish_color);
    txpos = xlims(1) + abs(xlims(2) - xlims(1))*0.1;
    typos = ylims(2) - abs(ylims(2) - ylims(1))*0.1;
    text(txpos,typos,sprintf('r=%.2f',c),'FontSize',rfont,'Color',reddish_color)
end