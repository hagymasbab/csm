function ge = gestaltGenerate(ge,N,varargin)
    parser = inputParser;
    addParamValue(parser,'precision',false,@islogical);
    addParamValue(parser,'verbose',true,@islogical);
    parse(parser,varargin{:});
    precision = parser.Results.precision; 
    verbose = parser.Results.verbose;
    if verbose
        fprintf('Generating synthetic data\n');
    end
    ge.N = N;
    if verbose
        fprintf('..Generating g values from a Dirichlet prior with concentration parameter %.2f\n',ge.sparsity);
    end
    ge.G = symmetricDirichlet(ge.sparsity,ge.k,ge.N);
    
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
        ge.X(n,:,:) = mvnrnd(means,ge.obsVar*eye(ge.Dx));
    end
    if verbose
        fprintf('\n');
    end
end
        