function patchDB = samplePatches(patchDim,trainingLength)

    % %%% Reset random seeds 
    % rand('state',sum(100*clock)); 
    % randn('state',sum(100*clock)); 
    numDims = patchDim^2;
    load(sprintf('wnodes_%d.mat',numDims));
    rawImages = load_images;
    % % making batches of image patches for training

    %trainingLength = 800*80;
    %batchLength = 800;

    fprintf('Generating batches of data...\n')

    %numBatches = trainingLength/batchLength;

    %batchData = zeros(numDims, batchLength, numBatches);
    %patchDB = zeros(numDims, trainingLength);
    %x = sample_images(rawImages, batchLength, patchDim);
    x = sample_images(rawImages, trainingLength, patchDim);
    x = log(x+2);
    % x = x - repmat(mean(x,1),batchLength,1);
    % x = reshape(x,batchLength,patchDim^2);
    x = x - repmat(mean(x,1),trainingLength,1);
    x = reshape(x,trainingLength,numDims);
    patchDB = do_whitening(x,wnode)';

    save(sprintf('bin/patches_vanhateren_%d.mat',numDims),'patchDB');
end
