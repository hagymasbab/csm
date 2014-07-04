function [cc,W,g_bias] = gestaltPretrain(ge,X,steps,varargin)
    % learn RBM weights between v and g as Gaussian units with CD1 
    
    parser = inputParser;
    addParamValue(parser,'alpha',0.01,@isnumeric);
    addParamValue(parser,'gstep',1,@isnumeric);
    addParamValue(parser,'cdstep',1,@isnumeric);
    addParamValue(parser,'initVar',0.1,@isnumeric);
    parse(parser,varargin{:});
    alpha = parser.Results.alpha;    
    gstep = parser.Results.gstep;    
    cdstep = parser.Results.cdstep;    
    initVar = parser.Results.initVar; 
    
    N = ge.N;
    % transform out batches
    if ndims(X) == 3
        X = reshape(X,size(X,1)*size(X,2),size(X,3));
        N = ge.N * ge.B;
    end
        
    % initialise W
    g_bias = randn(ge.k,1) * initVar;
    W = randn(ge.Dv,ge.k) * initVar;
    pA = pinv(ge.A);
%     data_corr = zeros(ge.Dv,ge.k);
%     fantasy_corr = zeros(ge.Dv,ge.k);
%     data_hiddenact = zeros(1,ge.k);
%     fantasy_hiddenact = zeros(1,ge.k);
%     samples = cell(1,steps);
    
    % transform each line of X into a V by the pseudoinverse of A
    V_data = pA * X'; % TODO check whether we have to transpose

    for s = 1:steps
        V = V_data;
        % take one sample for each v from fake-G
        G = gibbsG(V,W,gstep,g_bias);
        % record hidden activity
        data_hiddenact = sum(G,2);
        % record positive phase correlations        
        data_corr = V * G';
        
        for cds=1:cdstep
            % update V 
            V = W * G + randn(ge.Dv,N);
            % update G
            G = gibbsG(V,W,gstep,g_bias);
        end
        % record negative phase correlations
        fantasy_corr = V * G';
        % record hidden activity
        fantasy_hiddenact = sum(G,2);
        
        % update W
        W = W + alpha * (data_corr - fantasy_corr) / N;
        % update bias
        g_bias = g_bias + alpha * (data_hiddenact - fantasy_hiddenact) / N; 
        
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

function G = gibbsG(V,W,gstep,g_bias)
    k = size(W,2);
    N = size(V,2);
    G_mean = W' * V + repmat(g_bias,1,N);
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