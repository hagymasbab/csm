function lifetimeSparsity(cc,filterset,N,samples,burnin,dataset)
    Dx = size(cc{1},1);
    k = size(cc,2);
    
    ge = gestaltCreate('temp','Dx',Dx,'k',k,'B',1,'N',N,'filters',filterset, ...
        'obsVar',0.1,'g_shape',2,'g_scale',2,'z_shape',2,'z_scale',2,'nullComponent',false, ...
        'generateComponents',false,'generateData',false);
    
    datafile = sprintf('patches_%s_%d.mat',dataset,Dx);
    load(datafile);
    X = reshape(patchDB(:,1:N)',N,1,Dx);
    timings = ones(1,N) * (samples+burnin);
    [vs,gs,zs] = gestaltScheduling(X,timings,ge,1,0,true,gibbs,true);
    
    vdata = squeeze(vs(1,1,:,burnin+1:end,1,:)); % N x samples x Dx
    gdata = squeeze(gs(1,1,:,burnin+1:end,:)); % N x samples x k
    
    vrates = squeeze(mean(vdata,2));
    grates = squeeze(mean(gdata,2));
    
    % TODO calculate excess kurtosis
    
end