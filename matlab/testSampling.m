function all_G = testSampling(ge,L,sampler,prior)
    right = 0;
    all_G = zeros(ge.N,L,ge.k);
    rightindices = [];
    fprintf('Datum %d/',ge.N);
    for n=1:ge.N
        printCounter(n);
        s = gestaltGibbs(ge,n,L,'gSampler',sampler,'priorG',prior);
        all_G(n,:,:) = s(:,1:ge.k);
        meanG = mean(s(:,1:ge.k));
        [~,sampind] = max(meanG);
        [~,realind] = max(ge.G(n,:));
        if sampind == realind
            right = right + 1;
            rightindices = [rightindices n];
        else
            %ge.G(n,:)
        end
    end
    fprintf('\n%d/%d\n',right,ge.N);
    rightindices
end