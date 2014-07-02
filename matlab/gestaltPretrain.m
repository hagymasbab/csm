function cc = gestaltPretrain(ge,X,steps,alpha)
    % learn RBM weights between v and g as Gaussian units with CD1 
    
    % transform out batches
    if ndims(X) == 3
        X = reshape(X,size(X,1),size(X,2)*size(X,3));
    end
    
    % initialise W
    W = randn(ge.Dv,ge.k);
    pA = pinv(A);
    data_corr = zeros(ge.Dv,ge.k);
    fantasy_corr = zeros(ge.Dv,ge.k);
    
    for s = 1:steps
        % transform each line of X into a V by the pseudoinverse of A
        V = pA * X'; % TODO check whether we have to transpose
        % take one sample for each v from fake-G
        g = mvnrnd(W * v,eye(ge.k));
        
        % record positive phase correlations
        % update V 
        % update G
        % record negative phase correlations

        for n = ge.N
            % take one sample from V conditioned on X
            % now I just transform X by the pseudoinverse of A
            v = pA * X(n,b,:)';

            % take one sample from fake-G
            g = mvnrnd(W * v,eye(ge.k));

            % record positive phase correlations
            % update V 
            % update G
            % record negative phase correlations
        end
        
        % update W
        
    end
    
    % construct covariance components from W
    
end