function mode = sampleMode(samples,resolution)
    [counts,centres] = hist(samples,resolution);
    [~,maxidx] = max(counts);
    mode = centres(maxidx);
end