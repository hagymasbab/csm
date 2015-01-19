function [gRF,seeds] = gestaltGReceptiveFields(ge,cc,sampleNum)
    ge.cc = cc;
    ge.obsVar = 0.001;
    gRF = cell(1,ge.k);
    seeds = cell(1,ge.k);
    z = 1;
    for i = 1:ge.k
        act_G = zeros(ge.k,1);
        act_G(i,1) = 10;
        allX = zeros(sampleNum*ge.B,ge.Dx);
        for s = 1:sampleNum
            X = gestaltAncestralSample(ge,act_G,z);
            allX((s-1)*ge.B+1:s*ge.B,:) = X;            
        end
        %recField = zeros(ge.Dx,1);
        corrmat = corr(allX);
        corrmat = corrmat .* (ones(ge.Dx)-eye(ge.Dx));
        recField = max(corrmat);
        gRF{i} = recField;
        seeds{i} = sum(corrmat);
    end
end