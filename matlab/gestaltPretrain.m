function [cc,W] = gestaltPretrain(ge,X,steps,alpha)
    % learn RBM weights between v and g as Gaussian units with CD1 
    
    gstep = 3;
    N = ge.N;
    % transform out batches
    if ndims(X) == 3
        X = reshape(X,size(X,1)*size(X,2),size(X,3));
        N = ge.N * ge.B;
    end
        
    % initialise W
    W = randn(ge.Dv,ge.k);
    pA = pinv(ge.A);
    data_corr = zeros(ge.Dv,ge.k);
    fantasy_corr = zeros(ge.Dv,ge.k);
    samples = cell(1,steps);
    
    for s = 1:steps
        % transform each line of X into a V by the pseudoinverse of A
        V = pA * X'; % TODO check whether we have to transpose
        % take one sample for each v from fake-G
        G = gibbsG(V,W,gstep);
        % record positive phase correlations        
        data_corr = data_corr + V * G';
        
        % update V 
        V = W * G + randn(ge.Dv,N);
        % update G
        G = gibbsG(V,W,gstep);
        % record negative phase correlations
        fantasy_corr = fantasy_corr + V * G';
        
        % update W
        W = alpha * (data_corr - fantasy_corr) / N;
        
    end
    
    % construct covariance components from W
    cc = cell(1,ge.k);
    for k = 1:ge.k
        actc = zeros(ge.Dv);
        for i = 1:ge.Dv
            for j = 1:ge.Dv
                actc(i,j) = W(i,k) * W(j,k);
            end
        end
        cc{k} = actc;
    end
    
end

function G = gibbsG(V,W,gstep)
    k = size(W,2);
    N = size(V,2);
    G_mean = W' * V;
    G = G_mean + randn(k,N);
    % updating all G-s alternatingly with negative weight between each
    % other
    for i=1:gstep
        for j = 1:k
            restG = G;
            restG(j,:) = [];
            G(j,:) = G_mean(j,:) - sum(restG,1) + randn(1,N);
        end
    end
end