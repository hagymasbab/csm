function grad = gestaltParamGrad(ge,samples,cholesky)
    L = size(samples,2);
    N = size(samples,1);
    grad = cell(1,ge.k);
    cc = cell(1,ge.k);
    for i=1:ge.k
        grad{i} = zeros(ge.Dv);
        cc{i} = cholesky{i}' * cholesky{i};
    end
    fprintf('Sample %d/', N*L);
    for n=1:N
        for l=1:L
            printCounter((n-1)*L+l);
            g = squeeze(samples(n,l,1:ge.k));
            V = reshape(samples(n,l,ge.k+1:ge.k+ge.Dv*ge.B),ge.B,ge.Dv);
            VV = zeros(ge.Dv);
            for b=1:ge.B
                VV = VV + V(b,:)' * V(b,:);
            end
            iCv = inv(componentSum(g,cc));
            matr = ge.B * iCv - iCv * VV * iCv;
            for i=1:ge.k
                grad{i} = grad{i} - (g(i,1) * matr * cholesky{i}) / L;
            end
        end
    end
    fprintf('\n');
end
                