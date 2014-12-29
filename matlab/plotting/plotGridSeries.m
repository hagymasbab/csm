function plotGridSeries(data,redlines,verstrings,horstrings,vertitle,hortitle)
    figure();
    vernum = size(data,1);
    hornum = size(data,2);
    slen = size(data,4);
    if size(verstrings,2) == 0 
        for i = 1:vernum
            verstrings{end+1} = sprintf('%d',i);
        end
    end
    if size(horstrings,2) == 0 
        for i = 1:hornum
            horstrings{end+1} = sprintf('%d',i);
        end
    end
    ymin = min(data(:));
    ymax = max(data(:));
    for m = 1:vernum
        for kk = 1:hornum
            actdata = squeeze(data(m,kk,:,:));
            subplot(hornum,vernum,(kk-1)*vernum+m);
            plot(actdata');
            hold on;
            xlim([1 slen]);
            ylim([ymin ymax]);
            plot(mean(actdata)','LineWidth',3);            
            title(sprintf('%s %s %s %s',vertitle,verstrings{m},hortitle,horstrings{kk}),'FontSize',16);
            for t=1:length(redlines)
                act_t = redlines(t);
                plot([act_t;act_t],ylim(),'k-','LineWidth',2);
            end
        end
    end
end