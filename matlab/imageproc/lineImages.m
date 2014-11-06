function [patchDB,templates] = lineImages(N,Dx,k)
    templates = cell(1,k);
    temp_bool = cell(1,k);
    patchDB = zeros(N,Dx);
    imdim = sqrt(Dx);
    for i=1:k
        % TODO sample a random integer parameter set for lines
        a = randi([-2 2]);
        if a==0
            a=1;
        end
        if a < 0
            b = randi([imdim-2 imdim])
        else
            b = randi([0 2]);
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
        temp_bool{i} = logical(temp_img);
    end
    
    %coeffs = symmetricDirichlet(0.1,k,N); 
    coeffs = ones(N,k);
    for i=1:N
       img = 0.1 * randn(imdim,imdim);       
       % the order of occlusion should be random
       order = chooseKfromN(k,k);
       for j=1:k
           % select the pixels belonging to the templates 
           temppix = img(temp_bool{order(j)});
           % reduce their dynamic range to the extent of their coefficients
           meanpix = mean(temppix(:));
           img(temp_bool{order(j)}) = img(temp_bool{order(j)}) + max(coeffs(i,order(j)),1) * (meanpix - img(temp_bool{order(j)}));           
       end
       patchDB(i,:) = reshape(img,1,Dx);
    end
end