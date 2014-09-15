function [X,images,whitened] = gestaltImageStimulus(number,B,obsVar)
    % read the specified images
    im = imread('../natural_images/nat1.png');
    imdim = size(im,1);
    Dx = imdim^2;
    images = zeros(number,imdim,imdim);
    for i=1:number
        im = imread(sprintf('../natural_images/nat%d.png',i));        
        % transform the image to a real-valued [-5:5] vector from 0-255 uint8
        images(i,:,:) = double(im) / (255/10) - 5;
    end
    
    % whiten images
    images = reshape(images,number,Dx);
    whitened = whitenImages(images);             
    
    % add noise to each instance in the batch    
    X = permute(repmat(whitened,1,1,B),[1,3,2]) + obsVar * rand(number,B,Dx);
end
    