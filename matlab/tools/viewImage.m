function viewImage(data,varargin)
    p = inputParser;
    addParameter(p,'magnif',true,@islogical);
    addParameter(p,'usemax',false,@islogical);
    addParameter(p,'positive',false,@islogical);
    parse(p,varargin{:});
    magn = p.Results.magnif;
    usemax = p.Results.usemax;
    
    if iscell(data)
        if length(data) == 1
            data = data{1};
        else
            viewImageSet(data,'max',usemax);
            return;
        end
    end
    
    range = [-2 2];
    if usemax
        maxval = max(max(abs(data)));
        if maxval == 0
            maxval = 1;
        end
        range = [-maxval maxval];
    end
    if p.Results.positive
        range(1) = 0;
    end
    
    if ndims(data) > 2
        data = squeeze(data);
    end
    if min(size(data)) == 1
        imdim = sqrt(max(size(data)));
        data = reshape(data,imdim,imdim);
    end
    if magn
        IM = floor(24000 / size(data,1));
    else
        IM = 1;
    end
    imshow(data,'InitialMAgnification',IM,'colormap',gray,'DisplayRange',range);
end