subplot(1,2,1)
plot((0:size(batch_like,2)-1)',batch_like','LineWidth',1)
hold on
plot((0:size(batch_like,2)-1)',mean(batch_like),'LineWidth',3)
xlabel('Gradient ascent step #','FontSize',16)
ylabel('Log-likelihood on training sets','FontSize',16)
%set(gca,'XTick',0:size(batch_like,2)-1)
set(gca,'FontSize',16)
xlim([0 size(batch_like,2)-1])

subplot(1,2,2)
plot((0:length(test_like)-1)',test_like,'LineWidth',2)
xlabel('Gradient ascent batch #','FontSize',16)
ylabel('Log-likelihood on test set','FontSize',16)
if length(test_like) < 5
    set(gca,'XTick',0:length(test_like)-1)
end
set(gca,'FontSize',16)
xlim([0 length(test_like)-1])