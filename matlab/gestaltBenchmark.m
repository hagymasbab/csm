function diffs = gestaltBenchmark(ge,nRun,maxStep,name,hyperparams)
    defaults.samples                = 20;
%     defaults.increaseLikelihood     = false;
%     defaults.fullLikelihood         = false;
%     defaults.likelihoodSamples      = 10;
%     defaults.rateMethod             = 'componentwise_goal';
    defaults.learningRate           = 0.1;
%     defaults.multistep              = false;
    defaults.initCond               = 'empty';
    defaults.priorG                 = 'gamma';
    defaults.method                 = 'block';
    defaults.dataSource             = 'synthetic';
    
    defaults.dataPoints             = 50;
    defaults.batchSize              = 10;
    defaults.obsVar                 = 0.01;
    defaults.sparsity               = 0.2;
    
    if isempty(hyperparams)
        hyperparams = {{}};
    end
    nParams = size(hyperparams,2);
    parametrisations = cell(1,nParams);
    for hp=1:nParams
        fprintf('Parametrisation %d/%d ',nParams,hp);
        actparam = defaults;
        actparam = updateStruct(actparam,hyperparams{hp});
        parametrisations{hp} = actparam;
        save(sprintf('%s_params.mat',name),'parametrisations');
        
        if ~strcmp(actparam.dataSource,'synthetic')           
            load(actparam.dataSource); 
            if ~exist('patchDB')
                throw(MException('Gestalt:Benchmark:EmptySource','The data source %s does not contain a variable called patchDB',actparam.dataSource));
            end
            if ge.B ~= 1
                ge.B = 1;
                fprintf('Batch size is set to 1 as sampling conditioned in non-synthetic data is not implemented for observation batches.\n');
            end
            ge.N = actparam.dataPoints;
            ge.X = reshape(patchDB(:,1:actparam.dataPoints)',actparam.dataPoints,1,ge.Dv); % B needs to be 1
            ge.obsVar = actparam.obsVar;
        end
        
        fprintf('Run %d/',nRun);
        for r=1:nRun
            printCounter(r);
            if strcmp(actparam.dataSource,'synthetic')
                ge = gestaltGenerate(ge,actparam.dataPoints,'verbose',false,'batchSize',actparam.batchSize,'obsVar',actparam.obsVar,'sparsity',actparam.sparsity);
            end
            
            gestaltParameterEstimation(ge,ge.X,actparam.samples,maxStep,'shuffle','plot',0,'verbose',1, ... 
                'learningRate',actparam.learningRate,'initCond',actparam.initCond,'priorG',actparam.priorG, ... 
                'method',actparam.method,'syntheticData',strcmp(actparam.dataSource,'synthetic'));
            
            copyfile('iter.mat',sprintf('%s_iter_param%d_run%d.mat',name,hp,r));
        end
        fprintf('\n');
    end
    %plotDifferences(name,nRun,nParams);
end

function us = updateStruct(structure,cellArray)
    names = fieldnames(structure);
    us = structure;
    for i=1:size(names,1)        
        newValue = nextValue(cellArray,names(i));
        if ~isempty(newValue)
            us.(char(names(i))) = newValue;
        end
    end
end

function val = nextValue(cellArray,valName)
    val = '';
    for i=1:size(cellArray,2)
        if strcmp(cellArray{i},valName)
            val = cellArray{i+1};
            break;
        end
    end
end

    
    