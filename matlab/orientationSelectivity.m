function orientationSelectivity(nTrials,loadSamples,randseed)

    setrandseed(randseed);

    contrasts = [0.05 0.2 0.8];
    orient_shift = 10 * pi/180;
    central_orient = 0;
    orients = [central_orient-2*orient_shift central_orient-orient_shift central_orient central_orient+orient_shift central_orient+2*orient_shift];
    
    % TODO choose RF and cell
    
    
    stimuli = cell(1,length(contrasts)*length(orients));
    for o=1:length(orients)
        % TODO create Gabor
        for z=1:length(contrasts)
            % TODO put Gabor on gray background
        end
    end
    
    % TODO create model
    
    % TODO run scheduling
    
    % TODO calculate means
    
    % TODO plot 
    
end