function plotPair(leftData,rightData,plotRows,plotColumns,leftPlotIndex,plotMean,titleStrings,redLines,labels)
    subplot(plotRows,plotColumns,leftPlotIndex);
    plot(leftData');
    hold on;
    xlim([1,size(leftData,2)]);
    yl1 = ylim();
    if plotMean        
        plot(mean(leftData)','LineWidth',3);
    end
    title(titleStrings{1},'FontSize',16);
    
    subplot(plotRows,plotColumns,leftPlotIndex+1);
    plot(rightData');
    hold on;
    xlim([1,size(rightData,2)]);
    yl2 = ylim();
    if plotMean        
        plot(mean(rightData)','LineWidth',3);
    end
    title(titleStrings{2},'FontSize',16);
    if ~isempty(labels)
        legend(labels);
    end
    
    subplot(plotRows,plotColumns,leftPlotIndex);
    ylim([min(yl1(1),yl2(1)) max(yl1(2),yl2(2))]);
    plot([redLines{1};redLines{1}],ylim(),'r-');
    plot([redLines{1}+1;redLines{1}+1],ylim(),'k--');
    plot([size(rightData,2)-redLines{2};size(rightData,2)-redLines{2}],ylim(),'r-');
    plot([size(rightData,2)-redLines{2}-1;size(rightData,2)-redLines{2}-1],ylim(),'k--');
    plot([1; size(rightData,2)],[0;0],'k--');
    subplot(plotRows,plotColumns,leftPlotIndex+1);
    ylim([min(yl1(1),yl2(1)) max(yl1(2),yl2(2))]);
    plot([redLines{1};redLines{1}],ylim(),'r-');
    plot([redLines{1}+1;redLines{1}+1],ylim(),'k--');
    plot([size(rightData,2)-redLines{2};size(rightData,2)-redLines{2}],ylim(),'r-');
    plot([size(rightData,2)-redLines{2}-1;size(rightData,2)-redLines{2}-1],ylim(),'k--');
    plot([1; size(rightData,2)],[0;0],'k--');
end