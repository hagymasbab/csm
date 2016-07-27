function xw = do_whitening(x,wnode)

[nImage,imageSize] = size(x);

xw=zeros(size(x));

filterSize = size(wnode,1);

if mod(imageSize,filterSize) ~= 0
    error('Image size has to be an integer multiple of the filter size');
end
nStep = floor(sqrt(imageSize) / sqrt(filterSize));

for i=1:nImage
    if filterSize == imageSize
        xw(i,:) = (wnode*x(i,:)')';
    else
        actImg = reshape(x(i,:),sqrt(imageSize),sqrt(imageSize));
        wImg = zeros(sqrt(imageSize));
        for sx = 1:nStep
            for sy = 1:nStep               
                actPatch = actImg((sx-1)*sqrt(filterSize)+1:sx*sqrt(filterSize),(sy-1)*sqrt(filterSize)+1:sy*sqrt(filterSize));
                wPatch = reshape(wnode * actPatch(:),sqrt(filterSize),sqrt(filterSize));
                wImg((sx-1)*sqrt(filterSize)+1:sx*sqrt(filterSize),(sy-1)*sqrt(filterSize)+1:sy*sqrt(filterSize)) = wPatch;
            end
        end
        xw(i,:) = wImg(:)';
    end
end
