function [patchDB,templates] = lineImages(N,Dx,k)
    templates = cell(1,k);
    patchDB = zeros(N,Dx);
    imdim = sqrt(Dx);
    shifting = imdim/2;
    for i=1:k
        % sample a random integer parameter set for lines
        a = randi([-1 1]);
        if a==0
            a=1;
        end
        if a < 0
            b = randi([imdim-shifting imdim]);
        else
            b = randi([0 shifting]);
        end
        % the image is imagined to be in the upper-right quadrant
        temp_img = zeros(imdim);
        prev_val=b;
        for j=1:imdim
            val = a*j+b;
            val = max(val,1);
            val = min(val,imdim);        
            if prev_val <= val
                temp_img(j,prev_val+1:val) = 1;            
            else
                temp_img(j,val:prev_val) = 1;            
            end
            prev_val = val;
        end
        templates{i} = temp_img;
    end
    
    coeffs = symmetricDirichlet(0.1,k,N); 
    %coeffs = ones(N,k);
    for i=1:N
        img = createImageStimulus(templates,coeffs(i,:)');
        patchDB(i,:) = reshape(img,1,Dx);
    end
end