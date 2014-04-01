function viewImage(data,varargin)
    p = inputParser;
    addParamValue(p,'magnif',true,@islogical);
    parse(p,varargin{:});
    magn = p.Results.magnif;
    
    if min(size(data)) == 1
        imdim = sqrt(max(size(data)));
        data = reshape(data,imdim,imdim);
    end
    if magn
        IM = floor(24000 / size(data,1));
    else
        IM = 1;
    end
    imshow(data,'InitialMAgnification',IM,'colormap',jet,'DisplayRange',[-1 1]);
end