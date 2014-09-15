function whitened = whitenImages(images)
    % whiten images: http://www.mathworks.com/matlabcentral/fileexchange/34471-data-matrix-whitening/content/whiten.m
    % each image should be a row vector
    A = images'*images;
    [V,D,~] = svd(A);
    whMat = sqrt(size(images,1)-1)*V*sqrtm(inv(D + eye(size(D))*0.0001))*V';
    whitened = images*whMat;
end