function patchDB = samplePatches(patchDim,trainingLength,directory,name)

    numDims = patchDim^2;
    
    % load previously generated whitening filters
    load(sprintf('wnodes_%d.mat',numDims));
    
    if strcmp(directory,'vanhateren')
        rawImages = load_images; % this only works for Van Hateren images
    else
        files = dir(directory);
        files = files(~[files.isdir]);
        picNum = size(files,1);
        first_image = imread(strcat(directory,'/',files(1).name));
        rawImages = zeros(picNum,size(first_image,1),size(first_image,2));
        for i=1:picNum
            rawImages(i,:,:) = imread(strcat(directory,'/',files(i).name));
        end
    end
        
    % taking appropriately sized patches from the images randomly
    x = sample_images(rawImages, trainingLength, patchDim);
    
    %taking the logarithm of the pixel values
    x = log(x+2);
    
    % mean-centering the images
    x = x - repmat(mean(x,1),trainingLength,1);
    
    x = reshape(x,trainingLength,numDims);
    
    % whiten images
    patchDB = do_whitening(x,wnode)';

    save(sprintf('bin/patches_%s_%d.mat',name,numDims),'patchDB');
end
