function plotShiftedSeries(v)
    if ndims(v) > 2
        v = squeeze(v);
        if ndims(v) > 2
            error('more than 2 dims in data');
        end
    end
    nSeries = size(v,1);
    % step should be related to variance of the series
    step = max(abs(v(:)));
    %step = 1;
    shiftmat = repmat((0:step:(nSeries-1)*step)',1,size(v,2));
    h = plot((v+shiftmat)');
    c = get(h,'Color');
    axis off;
    hold on;
    xlim([1 size(v,2)]);
    means = mean(v,2);
    stds = std(v,0,2);
    for i = 1:length(means)
        out = -10;
        val = means(i) + (i-1)*step;
        mh = plot([out max(xlim)],[val val],'Color',c{i},'LineWidth',2);
        dh = plot([out out],[val-stds(i) val+stds(i)],'Color',c{i},'LineWidth',5);
        dh1 = plot([out out+2],[val-stds(i) val-stds(i)],'Color',c{i},'LineWidth',5);
        dh2 = plot([out out+2],[val+stds(i) val+stds(i)],'Color',c{i},'LineWidth',5);
        set(mh,'Clipping','off');
        set(dh,'Clipping','off');
        set(dh1,'Clipping','off');
        set(dh2,'Clipping','off');
    end
end