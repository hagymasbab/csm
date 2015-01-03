function samples = mergeSamples(vsamp,gsamp,zsamp)
    nSamp = size(vsamp,1);
    B = size(vsamp,2);
    Dv = size(vsamp,3);
    samples = [gsamp reshape(vsamp,[nSamp B*Dv]) zsamp];
end