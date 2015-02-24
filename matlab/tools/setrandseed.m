function setrandseed(randseed)
    if strcmp(randseed,'last')
        % TODO check if exsits
        load lastrandseed;
    elseif strcmp(randseed,'leave')
        return
    end
    s = RandStream('mt19937ar','Seed',randseed);
    RandStream.setGlobalStream(s);
    randseed = s.Seed;
    save('bin/lastrandseed.mat','randseed');
end