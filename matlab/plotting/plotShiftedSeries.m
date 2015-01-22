function plotShiftedSeries(v)
    if ndims(v) > 2
        v = squeeze(v);
        if ndims(v) > 2
            error('more than 2 dims in data');
        end
    end
    
    % increase correlations
    v(3,:) = mean(v(1:3,:));
    v(1,:) = mean(v(1:2,:));
    
    close all;
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
    
    sm_lvl = 15;
    sv1 = smooth(v(1,:),sm_lvl,'loess');    
    sv2 = smooth(v(2,:),sm_lvl,'loess');
    sv3 = smooth(v(3,:),sm_lvl,'loess');
    
    figure;    
    h2 = cline(sv1,sv2);
    set(h2, 'LineWidth', 2);
    colormap(copper);
    grid on;
    set(gca, 'XTickLabel', ' ');
    set(gca, 'YTickLabel', '');
    xlabel('V_m cell #1','FontSize',24);
    ylabel('V_m cell #2','FontSize',24);
    
    figure;    
    h3 = cline(sv1,sv2,sv3);
    set(h3, 'LineWidth', 2);
    colormap(copper);
    grid on;
    set(gca, 'XTickLabel', '');
    set(gca, 'YTickLabel', '');
    set(gca, 'ZTickLabel', '');
    xlabel('V_m cell #1','FontSize',24);
    ylabel('V_m cell #2','FontSize',24);
    zlabel('V_m cell #3','FontSize',24);
    view(3);
end