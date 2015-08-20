function filterPairCorrelations(corr1,orient1,corr2,orient2)
    nFilt = size(corr1,1);
    maxOr = 180;
    
    toPlot1 = [];
    toPlot2 = [];
    binidx1 = toBinindex(orient1,1:maxOr-1);
    binidx2 = toBinindex(orient2,1:maxOr-1);
    binned_average = nan(maxOr,maxOr);
    binned_variance = nan(maxOr,maxOr);
    
    binned_lists = cell(maxOr,maxOr);
    
    for i=1:nFilt
        for j=i+1:nFilt
           larger_c = max([orient1(i) orient1(j)]);
           smaller_c = min([orient1(i) orient1(j)]);
           toPlot1 = [toPlot1; larger_c smaller_c corr1(i,j)];
           larger_c = max([orient2(i) orient2(j)]);
           smaller_c = min([orient2(i) orient2(j)]);
           toPlot2 = [toPlot2; smaller_c larger_c corr2(i,j)];
           
           larger_bin = max([binidx1(i) binidx1(j)]);
           smaller_bin = min([binidx1(i) binidx1(j)]);
           binned_lists{larger_bin,smaller_bin} = [binned_lists{larger_bin,smaller_bin} corr1(i,j)];           
           larger_bin = max([binidx2(i) binidx2(j)]);
           smaller_bin = min([binidx2(i) binidx2(j)]);
           binned_lists{smaller_bin,larger_bin} = [binned_lists{smaller_bin,larger_bin} corr2(i,j)];
           
        end
    end
    
    for i=1:maxOr
        for j=1:maxOr
            binned_average(i,j) = mean(binned_lists{i,j});
            binned_variance(i,j) = var(binned_lists{i,j});
        end
    end
    
    close all
    
%     subplot(1,2,1);
    plotWithNan(binned_average)
    hold on
    plot([0 maxOr],[0 maxOr],'k','LineWidth',2);
    hold off
%     subplot(1,2,2);
%     plotWithNan(binned_variance);
%     hold on
%     plot([0 maxOr],[0 maxOr],'k','LineWidth',2);
%     hold off
    
    figure
    diff = abs(binned_average-binned_average');
    plotWithNan(diff./(abs(binned_average) + abs(binned_average')));
    
    figure
    scatter(abs(binned_average(:)),diff(:))
    
%     figure
%     scatter(abs(binned_average(:)),binned_variance(:));
%     
%     [~,indices] = sort(abs(toPlot1(:,3)));
%     toPlot1 = toPlot1(indices,:);
%     
%     [~,indices] = sort(abs(toPlot2(:,3)));
%     toPlot2 = toPlot2(indices,:);
   

%     figure 
%     scatter(toPlot1(:,1),toPlot1(:,2),10,toPlot1(:,3),'filled');
%     hold on;
%     scatter(toPlot2(:,1),toPlot2(:,2),10,toPlot2(:,3),'filled');
%     plot([0 maxOr],[0 maxOr],'k','LineWidth',2);    
%     xlim([0 maxOr]);
%     ylim([0 maxOr]);
%     hold off;
end

function binindices = toBinindex(v,lowerbounds)
    % first and last bins are assumed to be infinite
    binindices = zeros(size(v));
    for i=1:length(v)
        binindices(i) = length(lowerbounds)+1;
        for j=1:length(lowerbounds)
            if v(i) <= lowerbounds(j)
                binindices(i) = j;
                break;
            end
        end
    end
end

function plotWithNan(data)
    pcolor([data nan(size(data,2),1); nan(1,size(data,1)+1)]);
    shading flat;
    set(gca, 'ydir', 'reverse');
end