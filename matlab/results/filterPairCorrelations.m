function filterPairCorrelations(corr1,orient1,corr2,orient2)
    nFilt = size(corr1,1);
    
    toPlot1 = [];
    toPlot2 = [];
    
    for i=1:nFilt
        for j=i+    1:nFilt
           larger_c = max([orient1(i) orient1(j)]);
           smaller_c = min([orient1(i) orient1(j)]);
           toPlot1 = [toPlot1; larger_c smaller_c corr1(i,j)];
           larger_c = max([orient2(i) orient2(j)]);
           smaller_c = min([orient2(i) orient2(j)]);
           toPlot2 = [toPlot2; smaller_c larger_c corr2(i,j)];
        end
    end
%     
%     [~,indices] = sort(abs(toPlot1(:,3)));
%     toPlot1 = toPlot1(indices,:);
%     
%     [~,indices] = sort(abs(toPlot2(:,3)));
%     toPlot2 = toPlot2(indices,:);
   
    scatter(toPlot1(:,1),toPlot1(:,2),100,toPlot1(:,3),'filled');
    hold on;
    scatter(toPlot2(:,1),toPlot2(:,2),100,toPlot2(:,3),'filled');
    plot([0 180],[0 180],'k','LineWidth',2);    
    xlim([0 180]);
    ylim([0 180]);
    hold off;
end