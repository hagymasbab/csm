function ge = gestaltGenerate(ge,N,precision)
    fprintf('Generating synthetic data\n');
    ge.N = N;
    fprintf('..Generating g values from a Dirichlet prior with concentration parameter %.2f\n',ge.sparsity);
    ge.G = symmetricDirichlet(ge.sparsity,ge.k,ge.N);
    
    fprintf('..Generating v and x values %d/',ge.N);
    ge.V = zeros(ge.N,ge.B,ge.Dv);
    ge.X = zeros(ge.N,ge.B,ge.Dx);
    for n=1:ge.N 
        printCounter(n);
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
    fprintf('\n');
end
        