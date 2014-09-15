function viewImageSet(images)
    subplot = @(m,n,p) subtightplot (m, n, p, [0.025 0.001], [0 0.025], [0 0.01]);
    % each image should be a row
    N = size(images,1);
    % maximum images to plot is 16 x 16
    maxim = 12^2;
    if N > maxim;
        images = images(1:maxim,:);
        N = maxim;
    end
    xdim = floor(sqrt(N));
    ydim = ceil(N/xdim);
    
    clf;
    for i=1:N
        subplot(xdim,ydim,i);
        viewImage(images(i,:),'useMax',true);
    end
end