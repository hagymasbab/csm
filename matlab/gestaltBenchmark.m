function diffs = gestaltBenchmark(ge,nRun,maxStep,name,hyperparams)
    defaults.samples                = 20;
%     defaults.increaseLikelihood     = false;
%     defaults.fullLikelihood         = false;
%     defaults.likelihoodSamples      = 10;
%     defaults.rateMethod             = 'componentwise_goal';
    defaults.learningRate           = 1;
%     defaults.multistep              = false;
    defaults.initCond               = 'empty';
    defaults.priorG                 = 'gamma';
    defaults.method                 = 'block';
    
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
        
        fprintf('Run %d/',nRun);
        for r=1:nRun
            printCounter(r);
            ge = gestaltGenerate(ge,actparam.dataPoints,'verbose',false,'batchSize',actparam.batchSize,'obsVar',actparam.obsVar,'sparsity',actparam.sparsity);
            gestaltParameterEstimation(ge,ge.X,actparam.samples,maxStep,'shuffle','plot',0,'verbose',1,'learningRate',actparam.learningRate,'initCond',actparam.initCond,'priorG',actparam.priorG,'method',actparam.method);
            copyfile('iter.mat',sprintf('%s_iter_param%d_run%d.mat',name,hp,r));
        end
        fprintf('\n');
    end
    plotDifferences(name,nRun,nParams);
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

    
    