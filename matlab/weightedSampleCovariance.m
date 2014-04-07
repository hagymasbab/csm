function gVV = weightedSampleCovariance(samples,k,B)
    N = size(samples,1);
    L = size(samples,2);
    sdim = size(samples,3);
    Dv = (sdim - k) / B;
    gVV = cell(1,k);
    for j=1:k
        gVV{j} = zeros(Dv);
    end
    for n=1:N
        for l=1:L
            g = reshape(samples(n,l,1:k),1,k);
            v_batch = reshape(samples(n,l,k+1:sdim),B,Dv); 
            vv = zeros(Dv);
            for b=1:B
                v = v_batch(b,:);
                vv = vv + v'*v;
            end
            vv = vv / B;
            for j=1:k
                gVV{j} = gVV{j} + (1/(L*N)) * g(1,j) * vv;
            end
        end
    end
end