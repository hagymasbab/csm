function patchDB = createPatchDB(path,patchSize,mat_file)
    % mat file is a boolean specifying whether the images are in a
    % preprocessed .mat file or some generic image format like .png
    files = dir(path);
    files = files(~[files.isdir]);
    picNum = size(files,1);
    patches = [];
    for i=1:picNum
        printCounter(i,'maxVal',picNum,'stringVal','Picture');
        % load image
        if mat_file
            load(strcat(path,'/',files(i).name));
            act_image = LUM_Image;  % this is of course highly specific to the vanhateren image database, change if needed
        else
            act_image = imread(strcat(path,'/',files(i).name));
        end        
        
        px = floor(size(act_image,1) / patchSize);
        py = floor(size(act_image,2) / patchSize);
                
        for x = 1:px
            for y = 1:py
                % select patch
                act_patch = act_image(1+(x-1)*patchSize:x*patchSize,1+(y-1)*patchSize:y*patchSize);
                patches = [patches; act_patch(:)'];
            end
        end
    end    
        
    % whiten patches 
    patchDB = whitenImages(patches);
    %patchDB = patches;
    save(sprintf('patches_%d.mat',patchSize^2),'patchDB');
    
end