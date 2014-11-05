function [patchDB,templates] = lineImages(N,Dx,k)
    templates = cell(1,k);
    temp_bool = cell(1,k);
    patchDB = zeros(N,Dx);
    imdim = sqrt(Dx);
    for i=1:k
        % TODO sample a random integer parameter set for lines
        a = randi();
        b = randi();
        % the image is imagined to be in the upper-right quadrant
        temp_img = zeros(imdim);
        for j=1:imdim
            val = a*j+b;
            if val<=imdim && val>0
                temp_img(j,val) = 1;
            end
        end
        templates(i) = temp_img;
        temp_bool{i} = logical(temp_img);
    end
    
    coeffs = symmetricDirichlet(0.1,k,N); 
    for i=1:N
       img = 0.1 * randn(imdim,imdim);       
       % the order of occlusion should be random
       order = chooseKfromN(k,k);
       for j=1:k
           % select the pixels belonging to the templates 
           temppix = img(temp_bool{order(j)});
           % reduce their dynamic range to the extent of their coefficients
           meanpix = mean(temppix(:));
           img(temp_bool{order(j)}) = img(temp_bool{order(j)}) + coeffs(i,order(j)) * (meanpix - img(temp_bool{order(j)}));           
       end
       patchDB(i,:) = reshape(img,1,Dx);
    end
end