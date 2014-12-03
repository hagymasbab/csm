function stimuli = ICstimuli()
    imdim = 32;
    Dx = imdim^2;    
    nOrient = 4;
    shift = nOrient / 2;
    nRF = Dx / nOrient;
    rfsinarow = imdim/shift;
    nContrast = 5;
    Z = linspace(0.1,3,nContrast);
    % create a set of 32x32 Gabor filters with 4 orientations (1-vertical, 2-slash, 3-horizontal, 4-backslash)
    A = gaborFilterBank(imdim,imdim,2,2,[0;pi/4;pi/2;3*pi/4],[4]);
    % choose one orientation for each location
    locOrients = randi(nOrient,[nRF,1]);
    
    backgroundZ = 0.5;
    % choose an IC location
    ic_idx = floor(nRF/2 + rfsinarow/2);
    % define gestalts as sets of RF locations and orientation indices
    % now they have to be of equal length
    sh1 = 3; % single shift 
    sh2 = 5; % double shift
    templates = [ic_idx-sh2                 ic_idx-sh1                 ic_idx+sh1                 ic_idx+sh2;                 ...
                 ic_idx-(sh2*rfsinarow)+sh2 ic_idx-(sh1*rfsinarow)+sh1 ic_idx+(sh1*rfsinarow)-sh1 ic_idx+(sh2*rfsinarow)-sh2; ...
                 ic_idx-(sh2*rfsinarow)     ic_idx-(sh1*rfsinarow)     ic_idx+(sh1*rfsinarow)     ic_idx+(sh2*rfsinarow);     ...
                 ic_idx-(sh2*rfsinarow)-sh2 ic_idx-(sh1*rfsinarow)-sh1 ic_idx+(sh1*rfsinarow)+sh1 ic_idx+(sh2*rfsinarow)+sh2];
        
    stimuli = cell(nContrast,nOrient);
    for zidx = 1:nContrast 
        for o=1:nOrient
            act_orients = locOrients;        
            act_zs = ones(nRF,1) * backgroundZ;
            for i=1:size(templates,2)
                act_orients(templates(o,i)) = o;
                act_zs(templates(o,i)) = Z(zidx);
            end

            % set coeffs
            coeffs = zeros(Dx,1);
            for i=1:nRF
                if i ~= ic_idx
                    coeffs((i-1)*nOrient + act_orients(i,1)) = act_zs(i);
                end
            end
            stimuli{zidx,o} = A * coeffs;
        end
    end
    viewImageSet(stimuli,'max',false);
end
    