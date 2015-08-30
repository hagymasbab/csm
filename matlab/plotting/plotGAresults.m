close all
numcol = 2;
if exist('trueCC') && ~isempty(trueCC)    
    numcol = 3;    
end
numbatch = length(test_like)-1;
numstep = size(batch_like,2)-1;

batch_like_adjusted = zeros(size(batch_like));
if isfield(state_sequence{1},'batchSize')
    for b = 1:numbatch
        batch_like_adjusted(b,:) = batch_like_adjusted(b,:) / state_sequence{b}.batchSize;
    end
else
    batch_like_adjusted = batch_like ./ size(batch_indices,2);
end

subplot(2,numcol,1)
plot((0:numstep)',batch_like_adjusted','LineWidth',1)
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
if exist('trueLL') && trueLL ~= 0
    hold on;
    plot([0;numbatch],[trueLL/size(test_indices,2);trueLL/size(test_indices,2)],'r','LineWidth',2);
    hold off;
end
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
          norms(i,other) = norm(act_cc{other});
          dist = covcompRootMeanSquare(act_cc(kk),act_cc(other),1);
          %dist = dist / (norms(i,kk)*norms(i,other));
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

if exist('trueCC') && ~isempty(trueCC)    
    hold on;
    truenorms = {};
    for kk=1:ge.k
        act_truenorm = norm(trueCC{kk});
        truenorms{end+1} = act_truenorm;
        plot([0;numbatch],[act_truenorm;act_truenorm],'r','LineWidth',1);
    end
    hold off;
end

subplot(2,numcol,numcol+2)
plot((0:numbatch)',distances);
xlim([0 numbatch]);
xlabel('Gradient ascent batch #','FontSize',16)
%ylabel({'RMS distance of component pairs','normalised by the product of norms'},'FontSize',16)
ylabel('RMS distance of component pairs','FontSize',16)
set(gca,'FontSize',16)

if exist('trueCC') && ~isempty(trueCC)    
    hold on;
    for kk = 1:ge.k
      for other = kk+1:ge.k
          dist = covcompRootMeanSquare(act_cc(kk),act_cc(other),1);
          dist = dist / (truenorms{kk}*truenorms{other});
          plot([0;numbatch],[dist;dist],'r','LineWidth',1);
      end
    end
    hold off;
end

% difference from truth
if exist('trueCC') && ~isempty(trueCC)
    vcov = cov(reshape(ge.V,ge.N,ge.Dv));
    vcovdist_sum_avg = zeros(numbatch+1,1);
    truedist_sum_avg = zeros(numbatch+1,1);
    truedist_avg = zeros(numbatch+1,1);
    truedist_max = zeros(numbatch+1,1);
    truedist_sum_max = zeros(numbatch+1,1);
    true_cv = componentSum(1,trueCC);
    sigma_est = zeros(numbatch+1,1);
    differences = [];
    for i=1:numbatch+1
        act_cv = componentSum(1,state_sequence{i}.estimated_components);
        [truedist_sum_avg(i),~,truedist_sum_max(i)] = covcompRootMeanSquare(true_cv,act_cv,1);
        if ge.k < 4
            [truedist_avg(i),~,truedist_max(i)] = covcompRootMeanSquare(trueCC,state_sequence{i}.estimated_components,[]);
        end
        vcovdist_sum_avg(i) = covcompRootMeanSquare(vcov,act_cv,1);
        sigma_est(i) = state_sequence{i}.sigma_x;
        actdiff = upperTriangleValues(act_cv) - upperTriangleValues(true_cv);
        differences = [differences actdiff];
    end
    
    subplot(2,numcol,3)
    %plot((0:numbatch)',differences');
    hold on;
    if ge.k < 4
        plot((0:numbatch)',[truedist_sum_avg vcovdist_sum_avg truedist_avg],'LineWidth',2);    
    else
        plot((0:numbatch)',[truedist_sum_avg vcovdist_sum_avg],'LineWidth',2);
    end
    hold off;
    legend({'sum','vcov','perm'},'FontSize',12)
    xlim([0 numbatch]);
    xlabel('Gradient ascent batch #','FontSize',16)
    ylabel('Mean RMS distance from truth','FontSize',16)
    set(gca,'FontSize',16)
    
    subplot(2,numcol,numcol + 3)
    plot((0:numbatch)',sigma_est,'LineWidth',2);
    hold on;
    xlim([0 numbatch]);
    ylim([0 1]);
    plot([0;numbatch],[trueSigma;trueSigma],'r','LineWidth',2);
    hold off;
    xlabel('Gradient ascent batch #','FontSize',16)
    ylabel('Estimated and true obs. noise','FontSize',16)
    set(gca,'FontSize',16)        
end
