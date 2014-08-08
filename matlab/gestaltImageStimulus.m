function X = gestaltImageStimulus(number,B,obsVar)
    % read the specified image
    im = imread(sprintf('../natural_images/nat%d.png',number));
    % transform the image to a real-valued [-5:5] vector from 0-255 uint8
    imvec = double(im(:)) / (255/10) - 5;
    % add noise to each instance in the batch
    Dx = size(imvec,1);
    X = repmat(imvec',B,1) + obsVar * rand(B,Dx);
end
    