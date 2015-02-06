function orientationSelectivity(nTrials,loadSamples,randseed)

    setrandseed(randseed);
    Dx = 64;
    imSize = sqrt(Dx);
    k = 2;
    nSamples = 100;
    burnin = 50;

    contrasts = [0.05 0.2 0.8];
    orient_shift = 10 * pi/180;
    central_orient = 0;
    orients = [central_orient-2*orient_shift central_orient-orient_shift central_orient central_orient+orient_shift central_orient+2*orient_shift];        
    
    px = 0.5;
    py = 0.5;
    
    % TODO calculate cell id    
    cell_idx = 4;
    
    stimuli = cell(1,length(contrasts)*length(orients));
    for o=1:length(orients)
        % create Gabor      
        [~,~,act_gabor] = gabor('theta',orients(o),'lambda',4,'width',imSize,'height',imSize,'px',px,'py',py,'Sigma',1.8);
        for z=1:length(contrasts)
            % put Gabor on gray background
            act_stimulus = zeros(imSize) + contrasts(z) * act_gabor;
            stimuli{1,length(contrasts)*(o-1)+z} = act_stimulus(:);
        end
    end
    
    viewImageSet(stimuli,'max',false);
    
    % set timings
    timings = (nSamples+burnin) * ones(1,length(contrasts)*length(orients));
    
    % create model
    ge = gestaltCreate('temp','Dx',Dx,'k',k,'B',1,'N',1,'filters','gabor_4or','obsVar',1,'nullComponent',false,'generateComponents',true,'generateData',false);
    
    % run scheduling
    %[vsamp,gsamp,zsamp] = gestaltScheduling(stimuli,timings,{ge},nTrials,obsNoise,true,'gibbs');
    
    % TODO calculate means
    
    % TODO plot 
    
end