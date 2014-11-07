function gRF = gestaltGReceptiveFields(ge,cc,sampleNum)
    ge.cc = cc;
    gRF = cell(1,ge.k);
    z = 1;
    for i = 1:ge.k
        act_G = zeros(ge.k,1);
        act_G(i,1) = 1;
        allX = zeros(sampleNum*ge.B,ge.Dx);
        for s = 1:sampleNum
            X = gestaltAncestralSample(ge,act_G,z,false);
            allX((s-1)*ge.B+1:s*ge.B,:) = X;            
        end
        %recField = zeros(ge.Dx,1);
        covmat = corr(allX);
        recField = max(covmat .* (ones(ge.Dx)-eye(ge.Dx)));
        gRF{i} = recField;
    end
end