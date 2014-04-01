function ge = gestaltGenerate(ge,N)
    fprintf('Generating synthetic data\n');
    ge.N = N;
    fprintf('..Generating g values from a Dirichlet prior with concentration parameter %.2f\n',ge.sparsity);
    ge.G = symmetricDirichlet(ge.sparsity,ge.k,ge.N);
    
    fprintf('..Generating hidden values\n');
    ge.V = zeros(ge.N,ge.B,ge.Dv);
    for n=1:ge.N 
        %printCounter(n);
        Cv = componentSum(ge.G(n,:)',ge.cc);
        ge.V(n,:,:) = mvnrnd(zeros(ge.B,ge.Dv),Cv);                   
    end
    fprintf('..Generating observed values\n');
    ge.X = zeros(ge.N,ge.B,ge.Dx);
    for n=1:ge.N
        means = reshape(ge.V(n,:,:),ge.B,ge.Dv);
        means = means * ge.A';
        ge.X(n,:,:) = mvnrnd(means,ge.obsVar*eye(ge.Dx));
    end
end
        