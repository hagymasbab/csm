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
    
    % whiten images: http://xcorr.net/2013/04/30/whiten-images-in-matlab/
%     mX = bsxfun(@minus,images,mean(images)); %remove mean
%     fX = fft(fft(mX,[],2),[],3); %fourier transform of the images
%     spectr = sqrt(mean(abs(fX).^2)); %Mean spectrum
%     whitened = ifft(ifft(bsxfun(@times,fX,1./spectr),[],2),[],3); %whitened X
%     whitened = reshape(whitened,number,Dx);
    
    % whiten images: http://www.mathworks.com/matlabcentral/fileexchange/34471-data-matrix-whitening/content/whiten.m
    images = reshape(images,number,Dx);
    A = images'*images;
    [V,D,~] = svd(A);
    whMat = sqrt(size(images,1)-1)*V*sqrtm(inv(D + eye(size(D))*0.0001))*V';
    whitened = images*whMat;             
    
    % add noise to each instance in the batch
    
    X = permute(repmat(whitened,1,1,B),[1,3,2]) + obsVar * rand(number,B,Dx);
end
    