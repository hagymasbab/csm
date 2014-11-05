function gRF = gestaltGReceptiveFields(ge,sampleNum)
    gRF = cell(1,ge.k);
    z = 1;
    for i = 1:ge.k
        act_G = zeros(ge.k,1);
        act_G(i,1) = 1;
        sumX = zeros(ge.Dx,1);
        for s = 1:sampleNum
            X = gestaltAncestralSample(ge,act_G,z,false);
            X = reshape(mean(X,1),ge.Dx,1);
            sumX = sumX + X;
        end        
        gRF{i} = sumX / sampleNum;
    end
end