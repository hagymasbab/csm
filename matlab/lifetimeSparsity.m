function lifetimeSparsity(cc,filterset,N,dataset)
    Dx = size(cc{1},1);
    k = size(cc,2);
    
    ge = gestaltCreate('temp','Dx',Dx,'k',k,'B',1,'N',N,'filters',filterset, ...
        'obsVar',0.1,'g_shape',2,'g_scale',2,'z_shape',2,'z_scale',2,'nullComponent',false, ...
        'generateComponents',false,'generateData',false);
    
    datafile = sprintf('patches_%s_%d.mat',dataset,Dx);
    load(datafile);
    X = reshape(patchDB(:,1:N)',N,1,Dx);
    
end