function img = createImageStimulus(templates,coefficients)
   imdim = size(templates{1},1);
   k = size(templates,2);
   
   img = 0.1 * randn(imdim,imdim) - 1;       
   % the order of occlusion should be random
   order = chooseKfromN(k,k);
   for j=1:k
       % select the pixels belonging to the templates 
       temp_bool = logical(templates{order(j)});
       temppix = img(temp_bool);
       % reduce their dynamic range to the extent of their coefficients
       meanpix = mean(temppix(:));
       act_coeff = min(coefficients(order(j)),1);
       img(temp_bool) = img(temp_bool) + act_coeff * (meanpix - img(temp_bool)) + act_coeff;           
   end
end