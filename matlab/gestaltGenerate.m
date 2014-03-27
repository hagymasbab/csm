function ge = gestaltGenerate(ge,N)
    fprintf('Generating synthetic data\n');
    ge.N = N;
    fprintf('..Generating g values from a Dirichlet prior with concentration parameter %.2f\n',ge.sparsity);
    ge.G = symmetricDirichlet(ge.sparsity,ge.k,ge.N);
    
    fprintf('..Generating hidden values %d/',N);
    ge.V = zeros(ge.N,ge.Dv);
    for n=1:ge.N 
        printCounter(n);
        Cv = componentSum(ge.G(n,:)',ge.cc);
        ge.V(n,:) = mvnrnd(zeros(1,ge.Dv),Cv);
    end
    fprintf('\n..Generating observed values\n');
    ge.X = mvnrnd(ge.V*ge.A',ge.obsVar*eye(ge.Dx));
end
        