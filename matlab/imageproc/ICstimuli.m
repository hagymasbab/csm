function stimuli = ICstimuli()
    imdim = 32;
    Dx = imdim^2;    
    nOrient = 4;
    nRF = Dx / nOrient;
    nContrast = 2;
    Z = linspace(0,10,nContrast);
    % create a set of 32x32 Gabor filters with 4 orientations
    A = gaborFilterBank(imdim,imdim,2,2,[0;pi/4;pi/2;3*pi/4],[4]);
    % choose one orientation for each location
    locOrients = randi(nOrient,[nRF,1]);
    coeffs = zeros(Dx,1);
    backgroundZ = 1;
    % choose an IC location
    ic_idx = floor(nRF/2);
    for i=1:nRF
        if i ~= ic_idx
            coeffs((i-1)*nOrient + locOrients(i,1)) = backgroundZ;
        end
    end
    baseImage = A*coeffs;    
    
    stimuli = cell(nContrast,nOrient);
    for o=1:nOrient
    end
end
    