function [gRF,seeds,all_corrs] = gestaltGReceptiveFields(ge,cc,sampleNum,plothist)
    close all;
    ge.cc = cc;
    ge.obsVar = 0.001;
    gRF = cell(1,ge.k);
    seeds = cell(1,ge.k);
    transformed_loc = cell(1,ge.k);
    z = 1;
    if plothist
        load filtermatching
        vertnum = 3;        
    end
    
    all_corrs = cell(1,ge.k);
    all_vars = zeros(ge.k,ge.Dv);
    for i = 1:ge.k
        printCounter(i,'maxVal',ge.k,'stringVal','Component');
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
    if plothist
        figure;viewImageSet(all_corrs);
        corr_mean = componentSum(1/ge.k,all_corrs);
        figure;viewImage(corr_mean);
        %sum(sum(corr_mean<0))
        corr_std = zeros(ge.Dv);
        for i = 1:ge.k       
            corr_std = corr_std + (all_corrs{i}-corr_mean).^2;
        end
        corr_std = sqrt((1/ge.k)*corr_std);
        figure;viewImage(corr_std);
        var_mean = mean(all_vars,1)';
        var_std = std(all_vars,0,1)';
        figure;
        for i = 1:ge.k       
            transformed_loc{i} = exp(seeds{i}/5);
            
            act_var = diag(cc{i});
            select_idx = act_var > var_mean + var_std;
            subplot(ge.k,vertnum,(i-1)*vertnum+1);
            hist(orients(select_idx));
            if i==1;title('Variance');end;
            xlim([0 180]);           
                                    
            select_idx = [];
            for j = 1:ge.Dv
                for k = j+1:ge.Dv
%                     all_corrs{i}(j,k)
%                     (corr_mean(j,k) + corr_std(j,k))
%                     corr_mean(j,k)
%                     corr_std(j,k)
%                     all_corrs{i}(j,k) > (corr_mean(j,k) + corr_std(j,k))
%                     pause
                    if all_corrs{i}(j,k) > 0 && all_corrs{i}(j,k) > (corr_mean(j,k) + 2*corr_std(j,k))
                    %if all_corrs{i}(j,k) > 0 && all_corrs{i}(j,k) > corr_mean(j,k)                          
                        if ~any(select_idx==j);select_idx = [select_idx;j];end;
                        if ~any(select_idx==k);select_idx = [select_idx;k];end;
                    end
                end
            end       
            size(select_idx)
            subplot(ge.k,vertnum,(i-1)*vertnum+2);
            hist(orients(select_idx));
            if i==1;title('Positive');end;
            xlim([0 180]);           
            
            select_idx = [];
            for j = 1:ge.Dv
                for k = j+1:ge.Dv
                    if  all_corrs{i}(j,k) < 0 && all_corrs{i}(j,k) < corr_mean(j,k) - 2*corr_std(j,k)
                        if ~any(select_idx==j);select_idx = [select_idx;j];end;
                        if ~any(select_idx==k);select_idx = [select_idx;k];end;
                    end
                end
            end                                    
            subplot(ge.k,vertnum,(i-1)*vertnum+3);
            hist(orients(select_idx));
            if i==1;title('Negative');end;
            xlim([0 180]);    
            
%             select_idx = seeds{i} > mean([mean(seeds{i});max(seeds{i})]);
%             %sum(select_idx)
%             %ymax = 10;
%             %subplot(vertnum,ge.k,i);
%             subplot(ge.k,vertnum,(i-1)*vertnum+1);
%             hist(orients(select_idx));
%             if i==1;title('Orientations');end;
%             xlim([0 180]);           
%             %ylim([0 ymax]);
%             %subplot(vertnum,ge.k,i+ge.k);  
%             subplot(ge.k,vertnum,(i-1)*vertnum+2);
%             hist(maxX(select_idx'));
%             if i==1;title('Vertical location');end;
%             xlim([1 sqrt(ge.Dx)]);
%             %ylim([0 ymax]);
%             %subplot(vertnum,ge.k,i+2*ge.k);
%             subplot(ge.k,vertnum,(i-1)*vertnum+3);
%             hist(maxY(select_idx'));
%             if i==1;title('Horizontal location');end;
%             xlim([1 sqrt(ge.Dx)]);
%             %ylim([0 ymax]);
        end

        figure;
        viewImageSet(transformed_loc);
    end
end