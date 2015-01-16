function gestaltLearnNatural(code,imdim,embatch,samplesize,k,learningRate)    
    if imdim == 16
        % read patchDB
        load('sejnowski_patches_16.mat');
        %patchDB = patchDB(:,1:100);
        % read filter set
        %filterfile = 'OF_256.mat';
        filterfile = 'gabor_4or_16.mat';
    else
        fprintf('Not implemented');
        return;
    end
    
    % create model    
    ge = gestaltCreate('temp','Dx',imdim^2,'k',k,'B',1,'N',embatch, ...
        'filters',filterfile,'obsVar',1,'g_shape',1,'g_scale',1,'z_shape',1,'z_scale',0.1,'nullComponent',false);        
        
    % set initial conditions
    if code == 0
        initCond = 'empty';
    else
        initCond = code;
    end        
    
    % start learning
    gestaltEM(ge,patchDB',embatch,10000,samplesize,'shuffle','syntheticData',false,'initCond',initCond,'learningRate',learningRate);
end