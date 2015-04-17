function [gRF,seeds,angstds,transformed_loc] = gestaltGReceptiveFields(ge,cc,sampleNum,plothist)
    close all;
    ge.cc = cc;
    ge.obsVar = 0.001;
    gRF = cell(1,ge.k);
    seeds = cell(1,ge.k);
    transformed_loc = cell(1,ge.k);
    z = 1;
    Dx = size(cc{1},1);
    load(sprintf('filtermatching_%d.mat',Dx));
    vertnum = 3; 
    
    all_corrs = cell(1,ge.k);
    all_vars = zeros(ge.k,ge.Dv);
    for i = 1:ge.k
        if plothist
            printCounter(i,'maxVal',ge.k,'stringVal','Component');
        end
        all_vars(i,:) = diag(cc{i})';
        all_corrs{i} = corrcov(cc{i});        
        act_G = zeros(ge.k,1);
        act_G(i,1) = 10;
        if ge.B > 1
            allX = zeros(sampleNum*ge.B,ge.Dx);
            for s = 1:sampleNum
                X = gestaltAncestralSample(ge,act_G,z);
                allX((s-1)*ge.B+1:s*ge.B,:) = X;            
            end
        else
            allX = gestaltAncestralSample(ge,act_G,z,'N',sampleNum);
        end
        %recField = zeros(ge.Dx,1);
        corrmat = corr(allX);
        corrmat = corrmat .* (ones(ge.Dx)-eye(ge.Dx));        
        recField = max(corrmat);
        gRF{i} = recField;
        seeds{i} = sum(abs(corrmat));
    end
        
    corr_mean = componentSum(1/ge.k,all_corrs);
    %figure;viewImage(corr_mean);
    %sum(sum(corr_mean<0))
    corr_std = zeros(ge.Dv);
    all_corrs = nodiag(all_corrs);
    for i = 1:ge.k       
        corr_std = corr_std + (all_corrs{i}-corr_mean).^2;
    end
    corr_std = sqrt((1/ge.k)*corr_std);
    %figure;viewImage(corr_std);
    var_mean = mean(all_vars,1)';
    var_std = std(all_vars,0,1)';
    
    plot_ang = min(ge.k,15);
    angstds = zeros(plot_ang,vertnum);
    for i = 1:ge.k       
        transformed_loc{i} = exp(seeds{i}/5);

        act_var = diag(cc{i});
        select_idx = act_var > var_mean + var_std;
        if sum(select_idx) == 0
            select_idx = (1:ge.Dv)';
        end
        angstds(i,1) = circ_rad2ang(circ_std(circ_ang2rad(orients(select_idx)))) ;
        
        if plothist && i <= plot_ang
            if i==1
                figure;viewImageSet(all_corrs);
                figure;
            end
            subplot(plot_ang,vertnum,(i-1)*vertnum+1);
            hist(orients(select_idx));
            if i==1;title('Variance');end;
            xlim([0 180]);                
            %xlabel(sprintf('Angular STD %f',angstds(i,1)));
            set(gca,'XTickLabel',{});
        end
        
%         select_idx = [];
%         for j = 1:ge.Dv
%             for k = j+1:ge.Dv
%                 if all_corrs{i}(j,k) > 0 && all_corrs{i}(j,k) > (corr_mean(j,k) + corr_std(j,k))
%                     if ~any(select_idx==j);select_idx = [select_idx;j];end;
%                     if ~any(select_idx==k);select_idx = [select_idx;k];end;
%                 end
%             end
%         end
%         if sum(select_idx) == 0
%             select_idx = (1:ge.Dv)';
%         end
        maxel = 50;
        [~,high_x,high_y] = maxNElements(all_corrs{i},maxel,[]);
        select_idx = [];
        for j = 1:length(high_x)
            if ~any(select_idx==high_x(j));select_idx = [select_idx;high_x(j)];end;
            if ~any(select_idx==high_y(j));select_idx = [select_idx;high_y(j)];end;
        end
        angstds(i,2) = circ_rad2ang(circ_std(circ_ang2rad(orients(select_idx))));
        
        if plothist && i <= plot_ang
            subplot(plot_ang,vertnum,(i-1)*vertnum+2);
            hist(orients(select_idx));
            if i==1;title('Positive');end;
            xlim([0 180]);                   
            %xlabel(sprintf('Angular STD %f',angstds(i,2)));
            set(gca,'XTickLabel',{});
        end

%         select_idx = [];
%         for j = 1:ge.Dv
%             for k = j+1:ge.Dv
%                 if  all_corrs{i}(j,k) < 0 && all_corrs{i}(j,k) < corr_mean(j,k) - corr_std(j,k)
%                     if ~any(select_idx==j);select_idx = [select_idx;j];end;
%                     if ~any(select_idx==k);select_idx = [select_idx;k];end;
%                 end
%             end
%         end    
%         if sum(select_idx) == 0
%             select_idx = (1:ge.Dv)';
%         end
        [~,high_x,high_y] = maxNElements(-all_corrs{i},maxel,[]);
        select_idx = [];
        for j = 1:length(high_x)
            if ~any(select_idx==high_x(j));select_idx = [select_idx;high_x(j)];end;
            if ~any(select_idx==high_y(j));select_idx = [select_idx;high_y(j)];end;
        end
        angstds(i,3) = circ_rad2ang(circ_std(circ_ang2rad(orients(select_idx))));

        if plothist && i <= plot_ang
            subplot(plot_ang,vertnum,(i-1)*vertnum+3);
            hist(orients(select_idx));
            if i==1;title('Negative');end;
            xlim([0 180]);    
            %xlabel(sprintf('Angular STD %f',angstds(i,3)));
            set(gca,'XTickLabel',{});
        end
    end
    if plothist
        figure;
        viewImageSet(transformed_loc);
    end
end