function patchDB = samplePatches(patchDim,trainingLength)

    numDims = patchDim^2;
    
    % load previously generated whitening filters
    load(sprintf('wnodes_%d.mat',numDims));
    
    rawImages = load_images; % this only works for Van Hateren images

    % taking appropriately sized patches from the images randomly
    x = sample_images(rawImages, trainingLength, patchDim);
    
    %taking the logarithm of the pixel values
    x = log(x+2);
    
    % mean-centering the images
    x = x - repmat(mean(x,1),trainingLength,1);
    
    x = reshape(x,trainingLength,numDims);
    
    % whiten images
    patchDB = do_whitening(x,wnode)';

    save(sprintf('bin/patches_vanhateren_%d.mat',numDims),'patchDB');
end
