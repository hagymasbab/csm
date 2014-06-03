function testSampling(ge,L)
    right = 0;
    for n=1:ge.N
        s = gestaltGibbs(ge,n,L);
        meanG = mean(s(:,1:2));
        [~,sampind] = max(meanG);
        [~,realind] = max(ge.G(n,:));
        if sampind == realind
            right = right + 1;
        else
            %ge.G(n,:)
        end
    end
    fprintf('%d/%d\n',right,ge.N);
end