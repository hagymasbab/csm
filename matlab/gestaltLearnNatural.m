function gestaltLearnNatural(code,dataset,filterset,Dx,embatch,samplesize,burnin,k,nullcomp,maxStep,skipCheck,likelihood,learningRate)    

    % TODO try to create the patch DB and filters if needed
    
    datafile = sprintf('patches_%s_%d.mat',dataset,Dx);
    load(datafile);
    
    % TODO subsample patchDB to leave out a test set for crossvalidation
    
%     filterfile = sprintf('filters_%s_%d.mat',filterset,Dx);

%     if imdim == 16
%         % read patchDB
%         load('sejnowski_patches_16.mat');
%         %patchDB = patchDB(:,1:100);
%         % read filter set
%         %filterfile = 'OF_256.mat';
%         if strcmp(filterset,'gabor')
%             filterfile = 'gabor_4or_16.mat';
%         elseif strcmp(filterset,'OF')
%             filterfile = 'OF_256.mat';
%         end
%     else
%         fprintf('Not implemented');
%         return;
%     end
    
    % create model    
    ge = gestaltCreate('temp','Dx',Dx,'k',k,'B',1,'N',embatch, ...
        'filters',filterset,'obsVar',1,'g_shape',2,'g_scale',2,'z_shape',2,'z_scale',2,'nullComponent',nullcomp,'generateComponents',false,'generateData',false);        
        
    % set initial conditions
    if code == 0
        initCond = 'empty';
    else
        initCond = code;
    end        
    
    % this means automatically adapting the LR after the first gradient computation 
    %learningRate = 0;
    
    % start learning
    gestaltEM(ge,patchDB',embatch,maxStep,samplesize,'shuffle','syntheticData',false,'initCond',initCond,'learningRate',learningRate,'burnin',burnin,'computeLikelihood',likelihood,'verbose',2,'skipCheck',skipCheck);
end