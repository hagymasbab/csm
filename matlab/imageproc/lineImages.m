function [patchDB,templates] = lineImages(N,Dx,k)
    templates = cell(1,k);
%     temp_bool = cell(1,k);
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
%         temp_bool{i} = logical(temp_img);
    end
    
    coeffs = symmetricDirichlet(0.1,k,N); 
    %coeffs = ones(N,k);
    for i=1:N
        img = createImageStimulus(templates,coeffs(i,:)');
%        img = 0.1 * randn(imdim,imdim) - 1;       
%        % the order of occlusion should be random
%        order = chooseKfromN(k,k);
%        for j=1:k
%            % select the pixels belonging to the templates 
%            temppix = img(temp_bool{order(j)});
%            % reduce their dynamic range to the extent of their coefficients
%            meanpix = mean(temppix(:));
%            act_coeff = min(coeffs(i,order(j)),1);
%            img(temp_bool{order(j)}) = img(temp_bool{order(j)}) + act_coeff * (meanpix - img(temp_bool{order(j)})) + act_coeff;           
%        end
       patchDB(i,:) = reshape(img,1,Dx);
    end
end