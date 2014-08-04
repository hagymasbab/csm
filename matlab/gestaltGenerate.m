function ge = gestaltGenerate(ge,N,varargin)
    parser = inputParser;
    addParamValue(parser,'precision',false,@islogical);
    addParamValue(parser,'verbose',true,@islogical);
    addParamValue(parser,'batchSize',10,@isnumeric);
    addParamValue(parser,'obsVar',0.1,@isnumeric);
    addParamValue(parser,'sparsity',0.2,@isnumeric);
    parse(parser,varargin{:});
    precision = parser.Results.precision; 
    verbose = parser.Results.verbose;
    ge.B = parser.Results.batchSize;
    ge.obsVar = parser.Results.obsVar;
    ge.sparsity = parser.Results.sparsity;
    ge.N = N;
    
    if verbose
        fprintf('Generating synthetic data\n');
        fprintf('..Generating g values from a Dirichlet prior with concentration parameter %.2f\n',ge.sparsity);
    end        
    ge.G = symmetricDirichlet(ge.sparsity,ge.k,ge.N);
    
    if ge.contrast
        % sample z from a Gamma prior
        ge.Z = gamrnd(ge.z_shape,ge.z_scale,ge.N,1);
    end
    
    if verbose
        fprintf('..Generating v and x values %d/',ge.N);
    end
    ge.V = zeros(ge.N,ge.B,ge.Dv);
    ge.X = zeros(ge.N,ge.B,ge.Dx);
    for n=1:ge.N 
        if verbose
            printCounter(n);
        end
        if ~precision
            Cv = componentSum(ge.G(n,:)',ge.cc);
        else
            Cv = inv(componentSum(ge.G(n,:)',ge.pc));
        end
        ge.V(n,:,:) = mvnrnd(zeros(ge.B,ge.Dv),Cv);                                
        
        
        means = reshape(ge.V(n,:,:),ge.B,ge.Dv);
        means = means * ge.A';
        if ge.contrast
            means = ge.Z(n,1) * means;
        end
        ge.X(n,:,:) = mvnrnd(means,ge.obsVar*eye(ge.Dx));
    end
    if verbose
        fprintf('\n');
    end
end
        