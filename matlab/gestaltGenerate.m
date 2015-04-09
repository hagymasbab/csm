function ge = gestaltGenerate(ge,N,varargin)
    parser = inputParser;
    addParameter(parser,'precision',false,@islogical);
    addParameter(parser,'verbose',false,@islogical);
    addParameter(parser,'batchSize',10,@isnumeric);
    addParameter(parser,'obsVar',0.1,@isnumeric);
    addParameter(parser,'sparsity',0.2,@isnumeric);
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
    
    % heuristically, we are better of with a dirichlet for small k
    if ge.k < 10
        %ge.G = symmetricDirichlet(ge.sparsity,ge.k,ge.N);
        ge.G = symmetricDirichlet(0.1,ge.k,ge.N);
    else        
        ge.G = gestaltSamplePriorG(ge,N);
%         ge.G = zeros(ge.N,ge.k);
%         for i = 1:N
%             ge.G(i,:) = gestaltSamplePriorG(ge,1)';
%         end
    end
    
    if ge.contrast
        % sample z from a Gamma prior
        ge.Z = gamrnd(ge.z_shape,ge.z_scale,ge.N,1);
    end
    
    ge.V = zeros(ge.N,ge.B,ge.Dv);
    ge.X = zeros(ge.N,ge.B,ge.Dx);
    for n=1:ge.N 
        if verbose
            printCounter(n,'stringVal','..Generating v and x values','maxVal',ge.N,'newLine',true);
        end
        
        if ge.contrast
            z = ge.Z(n,1);
        else
            z = 1;
        end
        
        [ge.X(n,:,:),ge.V(n,:,:)] = gestaltAncestralSample(ge,ge.G(n,:)',z);
    end
end
        