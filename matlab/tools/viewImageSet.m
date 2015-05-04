function viewImageSet(images,varargin)
    p = inputParser;
    addParameter(p,'sum',false,@islogical);
    addParameter(p,'max',true,@islogical);  
    addParameter(p,'setmax',2,@isnumeric);  
    addParameter(p,'positive',false,@islogical);
    addParameter(p,'titles',{});
    parse(p,varargin{:});
    
    subplot = @(m,n,p) subtightplot (m, n, p, [0.025 0.001], [0 0.025], [0 0.01]);
    forcexdim = 0;
    forceydim = 0;
    if iscell(images)
        temp = [];
        % we assume it to be a row cell
        if size(images,1) > 1 && size(images,2) > 1
            forcexdim = size(images,1);
            forceydim = size(images,2);
        end
        for r=1:size(images,1)
            for c=1:size(images,2)
                temp = [temp; images{r,c}(:)'];
            end
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
    if forcexdim > 0
        xdim = forcexdim;
    end
    ydim = ceil(N/xdim);
    
    clf;
    for i=1:N
        subplot(xdim,ydim,i);
        viewImage(images(i,:),'useMax',p.Results.max,'positive',p.Results.positive,'setmax',p.Results.setmax);
        if size(p.Results.titles,2) >= i
            title(p.Results.titles{i});
        end
    end
end