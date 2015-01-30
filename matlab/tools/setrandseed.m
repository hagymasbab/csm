function setrandseed(randseed)
    if strcmp(randseed,'last')
        % TODO check if exsits
        load lastrandseed;
    end
    s = RandStream('mt19937ar','Seed',randseed);
    RandStream.setGlobalStream(s);
    randseed = s.Seed;
    save('lastrandseed.mat','randseed');
end