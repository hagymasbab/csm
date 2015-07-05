numcol = 2;
if exist('trueCC') && ~isempty(trueCC)
    numcol = 3;
end
numbatch = length(test_like)-1;
numstep = size(batch_like,2)-1;

subplot(2,numcol,1)
plot((0:numstep)',batch_like'/size(batch_indices,2),'LineWidth',1)
%hold on
%plot((0:size(batch_like,2)-1)',mean(batch_like),'LineWidth',3)
xlabel('Gradient ascent step #','FontSize',16)
ylabel('Log-likelihood per image on training sets','FontSize',16)
%set(gca,'XTick',0:size(batch_like,2)-1)
set(gca,'FontSize',16)
xlim([0 numstep])
yl1 = ylim();

subplot(2,numcol,2)
plot((0:numbatch)',test_like/size(test_indices,2),'LineWidth',2)
xlabel('Gradient ascent batch #','FontSize',16)
ylabel('Log-likelihood per image on test set','FontSize',16)
if length(test_like) < 5
    set(gca,'XTick',0:numbatch)
end
set(gca,'FontSize',16)
xlim([0 numbatch])
yl2 = ylim();
ylim([min([yl1(1) yl2(1)]) max([yl1(2) yl2(2)])]);
subplot(2,numcol,1)
ylim([min([yl1(1) yl2(1)]) max([yl1(2) yl2(2)])]);

% norms and distances
norms = zeros(numbatch+1,ge.k);
distances = zeros(numbatch+1,ge.k*(ge.k-1)/2);
for i=1:numbatch+1
    act_cc = state_sequence{i}.estimated_components;
    met_idx = 1;
    for kk = 1:ge.k
      norms(i,kk) = norm(act_cc{kk});
      for other = kk+1:ge.k
          dist = covcompRootMeanSquare(act_cc(kk),act_cc(other),1);
          distances(i,met_idx) = dist;
          met_idx = met_idx+1;
      end
    end
end

subplot(2,numcol,numcol+1)
plot((0:numbatch),norms);
xlim([0 numbatch]);
xlabel('Gradient ascent batch #','FontSize',16)
ylabel('2-norm of covariance components','FontSize',16)
set(gca,'FontSize',16)

subplot(2,numcol,numcol+2)
plot((0:numbatch)',distances);
xlim([0 numbatch]);
xlabel('Gradient ascent batch #','FontSize',16)
ylabel('RMS distance of component pairs','FontSize',16)
set(gca,'FontSize',16)

% difference from truth
if numcol == 3
    truedist_sum = zeros(numbatch+1,1);
    true_cv = componentSum(1,trueCC);
    for i=1:numbatch+1
        act_cv = componentSum(1,state_sequence{i}.estimated_components);
        truedist_sum(i) = covcompRootMeanSquare(true_cv,act_cv,1);
    end
    
    subplot(2,numcol,3)
    plot((0:numbatch)',truedist_sum,'LineWidth',2);
    xlim([0 numbatch]);
    xlabel('Gradient ascent batch #','FontSize',16)
    ylabel('Summed RMS distance from truth','FontSize',16)
    set(gca,'FontSize',16)
end
