function [cc,winc,gbiasinc,vbiasinc] = gestaltPretrain(ge,steps,randseed,varargin)
    % learn RBM weights between v and g as Gaussian units with CD1 
    
    parser = inputParser;
    addParamValue(parser,'alpha',0.01,@isnumeric);
    addParamValue(parser,'gstep',0,@isnumeric);
    addParamValue(parser,'cdstep',1,@isnumeric);
    addParamValue(parser,'initVar',0.1,@isnumeric);
    addParamValue(parser,'gbias',true,@islogical);
    addParamValue(parser,'vbias',true,@islogical);
    addParamValue(parser,'plot',false,@islogical);
    addParamValue(parser,'verbose',false,@islogical);
    addParamValue(parser,'GDistrib','normal');
    addParamValue(parser,'saveCC',0,@isnumeric);
    parse(parser,varargin{:});
    params = parser.Results;    
    
    if strcmp(randseed,'last')
        load lastrandseed;
    end
    s = RandStream('mt19937ar','Seed',randseed);
    RandStream.setGlobalStream(s);
    randseed = s.Seed;
    save('lastrandseed.mat','randseed');
    
    X = ge.X;
    % transform out batches
    if ndims(X) == 3
        X = reshape(X,size(X,1)*size(X,2),size(X,3));
        N = ge.N * ge.B;
    end
    N = size(X,1);
        
    % initialise W
    g_bias = randn(ge.k,1) * params.initVar;
    v_bias = randn(ge.Dv,1) * params.initVar;
    W = randn(ge.Dv,ge.k) * params.initVar;
    
    if params.saveCC > 0
        savenum = floor(steps / params.saveCC) + 1;
        cc_iter = cell(savenum,ge.k);
        cc_iter(1,:) = w2cc(ge,W);        
        savecount = 2;
    end
    
%     data_corr = zeros(ge.Dv,ge.k);
%     fantasy_corr = zeros(ge.Dv,ge.k);
%     data_hiddenact = zeros(1,ge.k);
%     fantasy_hiddenact = zeros(1,ge.k);
%     samples = cell(1,steps);
    winc = zeros(ge.Dv * ge.k,steps);
    params.gbiasinc = zeros(ge.k,steps);
    params.vbiasinc = zeros(ge.Dv,steps);
    
    % transform each line of X into the V space
    %V_data = pinv(ge.A) * X';
    V_data = ge.A' * X';
    
    for s = 1:steps
        if params.verbose
            printCounter(s,'maxVal',steps,'stringVal','Step');
        end
        V = V_data;
        % take one sample for each v from fake-G
        G = gibbsG(V,W,params.gstep,params.gbias,g_bias,params.GDistrib);
        % record G activity
        data_gact = sum(G,2);
        % record V activity
        data_vact = sum(V,2);
        % record positive phase correlations        
        data_corr = correlate(V,G);
        
        for cds=1:params.cdstep
            % update V 
            if params.vbias
                V_bias_input = repmat(v_bias,1,N);
            else
                V_bias_input = 0;
            end
            V = W * G + V_bias_input + randn(ge.Dv,N);
            % update G
            G = gibbsG(V,W,params.gstep,params.gbias,g_bias,params.GDistrib);
        end
        % record negative phase correlations
        fantasy_corr = correlate(V,G);
        % record G activity
        fantasy_gact = sum(G,2);
        % record V activity
        fantasy_vact = sum(V,2);
        
        % update W      
        W = W + params.alpha * (data_corr - fantasy_corr);
        winc(:,s) = W(:)';
        % update bias
        if params.gbias
            g_bias = g_bias + params.alpha * (data_gact - fantasy_gact) / N; 
            params.gbiasinc(:,s) = g_bias;
        end
        if params.vbias
            v_bias = v_bias + params.alpha * (data_vact - fantasy_vact) / N; 
            params.vbiasinc(:,s) = v_bias;
        end
               
        if params.saveCC > 0 && rem(s,params.saveCC) == 0
            cc_iter(savecount,:) = w2cc(ge,W);        
            savecount = savecount + 1;
            save('bin/cc_pret_iter.mat','cc_iter','ge');
        end
        
    end
    
    
    % construct covariance components from W
    
    cc = w2cc(ge,W);
    
    % plot 
    if params.plot
        close all;
%     cc = cell(1,ge.k);
        for k = 1:ge.k
%         actc = zeros(ge.Dv);
%         for i = 1:ge.Dv
%             for j = 1:ge.Dv
%                 actc(i,j) = W(i,k) * W(j,k);
%             end
%         end
%         cc{k} = actc;
        
            figure;
            viewImage(cc{k},'usemax',true);
        end
    end        
end

function cc = w2cc(ge,W)    
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

function G = gibbsG(V,W,gstep,gbias,g_bias,distrib)
    k = size(W,2);
    N = size(V,2);
    if gbias
        G_bias_input = repmat(g_bias,1,N);
    else
        G_bias_input = 0;
    end
    G_mean = W' * V + G_bias_input;
    G = updateG(G_mean,N,distrib);
    %G = G_mean + randn(k,N);
    % updating all G-s alternatingly with negative weight between each
    % other
    for i=1:gstep
        for j = 1:k
            restG = G;
            restG(j,:) = [];
            %G(j,:) = G_mean(j,:) - sum(restG,1) + randn(1,N);
            G = updateG(G_mean(j,:) - sum(restG,1),N,distrib);
        end
    end
end

function G = updateG(mu,N,distrib)
    k = size(mu,1);
    if strcmp(distrib,'normal')
        G = mu + randn(k,N);
    elseif strcmp(distrib,'lognormal')
        G = lognrnd(mu,1,k,N);
    else
        exit(1);
    end
end

function corr = correlate(A,B) 
    N = size(A,2);
    % standardise
    A = A ./ repmat(std(A,0,2),1,N);
    B = B ./ repmat(std(B,0,2),1,N);
    % covariance
    corr = A * B';
    % normalise
    corr = corr / N;
end