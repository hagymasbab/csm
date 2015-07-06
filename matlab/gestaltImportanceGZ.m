function [gsamp,zsamp,weights] = gestaltImportanceGZ(x,ge,L,randseed)
    setrandseed(randseed);
    
    priorG = gestaltSamplePriorG(ge,L,'checkValues',false);
    priorZ = gamrnd(ge.z_shape,ge.z_scale,[L 1]);
    
    for i=1:L
        
    end
end