function gestaltLearnNatural(Dx,k,embatch,samplesize,varargin)    
    parser = inputParser;
    addParameter(parser,'code',0,@isnumeric);    
    addParameter(parser,'learningRate',0,@isnumeric); 
    addParameter(parser,'dataset','vanhateren'); 
    addParameter(parser,'filterset','OF'); 
    addParameter(parser,'nullcomp',false,@islogical);
    addParameter(parser,'skipCheck',false,@islogical);
    addParameter(parser,'datasize',0,@isnumeric);
    addParameter(parser,'datashift',0,@isnumeric);
    addParameter(parser,'burnin',50,@isnumeric);
    addParameter(parser,'likelihood','batch'); 
    addParameter(parser,'maxStep',10000,@isnumeric);
    addParameter(parser,'randseed','shuffle'); 
    parse(parser,varargin{:});        
    params = parser.Results;      
        
    % TODO try to create the patch DB and filters if needed
    
    datafile = sprintf('patches_%s_%d.mat',params.dataset,Dx);
    load(datafile);
    if params.datasize > 0
        patchDB = patchDB(:,1+params.datashift:params.datasize+params.datashift);
    end
    
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
        'filters',params.filterset,'obsVar',1,'g_shape',2,'g_scale',2,'z_shape',2,'z_scale',2,'nullComponent',params.nullcomp,'generateComponents',false,'generateData',false);        
        
    % set initial conditions
    if params.code == 0
        initCond = 'empty';
    else
        initCond = params.code;
    end        
    
    % this means automatically adapting the LR after the first gradient computation 
    %learningRate = 0;
    
    % start learning
    gestaltEM(ge,patchDB',embatch,params.maxStep,samplesize,params.randseed,'syntheticData',false,'initCond',initCond, ...
        'learningRate',params.learningRate,'burnin',params.burnin,'computeLikelihood',params.likelihood,'verbose',2,'skipCheck',params.skipCheck);
end