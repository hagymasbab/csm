function smallImages = sample_images(imgs, sampleNo, varargin)

% smallImages = sample_images(sampleNo, varargin)

args=varargin;
nargs = length(args);
if (nargs>0)
    smallImageSize=args{1};
else
    smallImageSize = 20;
end

xEdgeLength = size(imgs,2);
yEdgeLength = size(imgs,3);

imgNo = size(imgs,1);

samples = sampleNo;

xNormEdgeLength = xEdgeLength - smallImageSize;
yNormEdgeLength = yEdgeLength - smallImageSize;

upperLefts = floor(rand(samples,2) .* repmat([xNormEdgeLength, yNormEdgeLength],samples, 1)) + 1;

smallImages = zeros(samples, smallImageSize, smallImageSize);
for s = 1:samples
    imgID=floor(rand(1)*imgNo)+1;
    r1 = upperLefts(s,1);
    c1 = upperLefts(s,2);
    smallImages(s, :, :) = ...
        imgs(imgID, r1:(r1+smallImageSize-1), c1:(c1+smallImageSize-1)); %, ...
%            smallImageSize*smallImageSize, 1);
end
