function [vsamp,gsamp,zsamp] = splitSamples(samples,k,B)
    nSamp = size(samples,1);
    Dv = (size(samples,2) - k - 1) / B;
    gsamp = samples(:,1:k);
    vsamp = reshape(samples(:,k+1:end-1),[nSamp B Dv]);
    zsamp = samples(:,end);
end