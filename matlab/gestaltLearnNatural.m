function gestaltLearnNatural(code,imdim,embatch,samplesize)    
    if imdim == 16
        % read patchDB
        load('sejnowski_patches_16.mat');
        % read filter set
        filterfile = 'OF_256.mat';
    else
        fprintf('Not implemented');
        return;
    end
    
    % create model
    % TODO find reasons for parameter values
    ge = gestaltCreate('temp','Dx',imdim^2,'k',imdim^2,'B',1,'N',embatch, ...
        'filters',filterfile,'obsVar',1,'g_scale',2,'z_shape',1,'z_scale',0.1,'nullComponent',true);        
        
    % set initial conditions
    if code == 0
        initCond = 'empty';
    else
        initCond = code;
    end        
    
    % start learning
    gestaltEM(ge,patchDB',embatch,10000,samplesize,'shuffle','syntheticData',false,'initCond',initCond);
end