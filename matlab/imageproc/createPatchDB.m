function patchDB = createPatchDB(path,patchSize,maxOrigSize)
    %path = '/home/banmi/imageDB/sejnowski_images';
    %path = '/home/banmi/botswana/keves';
    files = dir(path);
    files = files(~[files.isdir]);
    %files = files(3:end);
    picNum = size(files,1);
    patches = [];
    fprintf('Picture %d/',picNum);
    for i=1:picNum
        printCounter(i);
        % load image
        mat_file = false;
        if mat_file
            load(strcat(path,'/',files(i).name));
            act_image = LUM_Image;
        else
            act_image = imread(strcat(path,'/',files(i).name));
        end

        % subsample image   
%         if maxOrigSize > 0
%             sub_step = ceil(max(size(act_image)) / maxOrigSize);
%         else
%             sub_step = 1;
%         end
%         subsampled = act_image(1:sub_step:end,1:sub_step:end);
        subsampled = act_image;
        
%         if i == 1
%             fprintf('Original dims %d %d, sampling step %d, new dims %d %d\n',size(LUM_Image,1),size(LUM_Image,2),sub_step,size(subsampled,1),size(subsampled,2));
%         end
    
        px = floor(size(subsampled,1) / patchSize);
        py = floor(size(subsampled,2) / patchSize);
                
        for x = 1:px
            for y = 1:py
                % select patch
                act_patch = subsampled(1+(x-1)*patchSize:x*patchSize,1+(y-1)*patchSize:y*patchSize);
                patches = [patches; act_patch(:)'];
            end
        end
    end    
    fprintf('\n');
        
    % whiten patches 
    patchDB = whitenImages(patches);
    %patchDB = patches;
    save(sprintf('sejnowski_patches_%d.mat',patchSize),'patchDB');
    
end