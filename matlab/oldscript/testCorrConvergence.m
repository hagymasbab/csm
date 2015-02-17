function testCorrConvergence(ge,cc,nTrial,nsamps)
    ge.cc = cc;
    ge.obsVar = 0.001;
    z = 1;
    allcorrs = cell(nTrial,length(nsamps),ge.k);
    for t=1:nTrial
        fprintf('Trial %d/%d ',nTrial,t);
        for s=1:length(nsamps)
            printCounter(s,'maxVal',length(nsamps),'stringVal','SamplingLength');
            for i = 1:ge.k                
                act_G = zeros(ge.k,1);
                act_G(i,1) = 10;
                allX = gestaltAncestralSample(ge,act_G,z,'N',nsamps(s));
                %recField = zeros(ge.Dx,1);
                corrmat = corr(allX);
                allcorrs{t,s,i} = corrmat;                
            end            
        end
        save('corrmats.mat','allcorrs','nsamps','ge','cc');
    end
end