function [gRF,seeds] = gestaltGReceptiveFields(ge,cc,sampleNum,plothist)
    %close all;
    ge.cc = cc;
    ge.obsVar = 0.001;
    gRF = cell(1,ge.k);
    seeds = cell(1,ge.k);
    transformed_loc = cell(1,ge.k);
    z = 1;
    if plothist
        load filtermatching
        vertnum = 3;
        figure;
    end
    for i = 1:ge.k
        printCounter(i,'maxVal',ge.k,'stringVal','Component');
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
        if plothist
            transformed_loc{i} = exp(seeds{i}/5);
            select_idx = seeds{i} > mean([mean(seeds{i});max(seeds{i})]);
            %sum(select_idx)
            %ymax = 10;
            %subplot(vertnum,ge.k,i);
            subplot(ge.k,vertnum,(i-1)*vertnum+1);
            hist(orients(select_idx));
            if i==1;title('Orientations');end;
            xlim([0 180]);           
            %ylim([0 ymax]);
            %subplot(vertnum,ge.k,i+ge.k);  
            subplot(ge.k,vertnum,(i-1)*vertnum+2);
            hist(maxX(select_idx'));
            if i==1;title('Vertical location');end;
            xlim([1 sqrt(ge.Dx)]);
            %ylim([0 ymax]);
            %subplot(vertnum,ge.k,i+2*ge.k);
            subplot(ge.k,vertnum,(i-1)*vertnum+3);
            hist(maxY(select_idx'));
            if i==1;title('Horizontal location');end;
            xlim([1 sqrt(ge.Dx)]);
            %ylim([0 ymax]);
        end
    end        
    if plothist
        figure;
        viewImageSet(transformed_loc);
    end
end