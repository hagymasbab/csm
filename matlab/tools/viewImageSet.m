function viewImageSet(images,varargin)
    p = inputParser;
    addParamValue(p,'sum',false,@islogical);
    addParamValue(p,'max',true,@islogical);
    parse(p,varargin{:});
    
    subplot = @(m,n,p) subtightplot (m, n, p, [0.025 0.001], [0 0.025], [0 0.01]);
    if iscell(images)
        temp = [];
        % we assume it to be a row cell
        for c=1:size(images,2)
            temp = [temp; images{c}(:)'];
        end
        images = temp;
    end
    
    % each image should be a row
    N = size(images,1);
    % maximum images to plot is 16 x 16
    maxim = 12^2;
    
    if p.Results.sum
        sumimage = sum(images);
        if N < maxim
            images = [images; sumimage];
            N = N + 1;
        else
            images(maxim,:) = sumimage;
        end
    end
    
    if N > maxim;
        images = images(1:maxim,:);
        N = maxim;
    end
    xdim = floor(sqrt(N));
    ydim = ceil(N/xdim);
    
    clf;
    for i=1:N
        subplot(xdim,ydim,i);
        viewImage(images(i,:),'useMax',p.Results.max);
    end
end