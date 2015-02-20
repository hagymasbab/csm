function [gRF,seeds,angstds,transformed_loc] = gestaltGReceptiveFields(ge,cc,sampleNum,plothist)
    close all;
    ge.cc = cc;
    ge.obsVar = 0.001;
    gRF = cell(1,ge.k);
    seeds = cell(1,ge.k);
    transformed_loc = cell(1,ge.k);
    z = 1;
    load filtermatching
    vertnum = 3; 
    
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
    
    %figure;viewImageSet(all_corrs);
    corr_mean = componentSum(1/ge.k,all_corrs);
    figure;viewImage(corr_mean);
    %sum(sum(corr_mean<0))
    corr_std = zeros(ge.Dv);
    for i = 1:ge.k       
        corr_std = corr_std + (all_corrs{i}-corr_mean).^2;
    end
    corr_std = sqrt((1/ge.k)*corr_std);
    %figure;viewImage(corr_std);
    var_mean = mean(all_vars,1)';
    var_std = std(all_vars,0,1)';

    angstds = zeros(ge.k,vertnum);
    for i = 1:ge.k       
        transformed_loc{i} = exp(seeds{i}/5);

        act_var = diag(cc{i});
        select_idx = act_var > var_mean + var_std;
        angstds(i,1) = circ_rad2ang(circ_std(circ_ang2rad(orients(select_idx)))) ;
        
        if plothist
            if i==1
                figure;
            end
            subplot(ge.k,vertnum,(i-1)*vertnum+1);
            hist(orients(select_idx));
            if i==1;title('Variance');end;
            xlim([0 180]);                
            xlabel(sprintf('Angular STD %f',angstds(i,1)));
        end

        select_idx = [];
        for j = 1:ge.Dv
            for k = j+1:ge.Dv
                if all_corrs{i}(j,k) > 0 && all_corrs{i}(j,k) > (corr_mean(j,k) + 2*corr_std(j,k))
                    if ~any(select_idx==j);select_idx = [select_idx;j];end;
                    if ~any(select_idx==k);select_idx = [select_idx;k];end;
                end
            end
        end       
        angstds(i,2) = circ_rad2ang(circ_std(circ_ang2rad(orients(select_idx)))) ;
        
        if plothist
            subplot(ge.k,vertnum,(i-1)*vertnum+2);
            hist(orients(select_idx));
            if i==1;title('Positive');end;
            xlim([0 180]);                   
            xlabel(sprintf('Angular STD %f',angstds(i,2)));
        end

        select_idx = [];
        for j = 1:ge.Dv
            for k = j+1:ge.Dv
                if  all_corrs{i}(j,k) < 0 && all_corrs{i}(j,k) < corr_mean(j,k) - 2*corr_std(j,k)
                    if ~any(select_idx==j);select_idx = [select_idx;j];end;
                    if ~any(select_idx==k);select_idx = [select_idx;k];end;
                end
            end
        end          
        angstds(i,3) = circ_rad2ang(circ_std(circ_ang2rad(orients(select_idx)))) ;

        if plothist
            subplot(ge.k,vertnum,(i-1)*vertnum+3);
            hist(orients(select_idx));
            if i==1;title('Negative');end;
            xlim([0 180]);    
            xlabel(sprintf('Angular STD %f',angstds(i,3)));
        end
    end
    if plothist
        figure;
        viewImageSet(transformed_loc);
    end
end