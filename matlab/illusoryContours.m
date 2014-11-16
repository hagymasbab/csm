function illusoryContours(nTrials,nSamples,try_k,k,Dx,filters)
    % reproducing the effect of responses to illusory and real contours 
    % from Lee & Nguyen, PNAS, 2001.
    close all;
    
    % model parameters for generation and sampling
    generating_sigma = 0.001;
    sampling_sigma = 1;
    g_scale = 2;
    z_shape = 1;
    z_scale = 0.1;
    v_sampler = 'direct';
    sample_z = true;
    
    % create model
    if strcmp(filters,'OF')
        filters = sprintf('OF_%d.mat',Dx);
    end
    ge = gestaltCreate('temp','Dx',Dx,'k',k,'B',1,'N',2, ...
        'filters',filters,'obsVar',generating_sigma,'g_scale',g_scale,'z_shape',z_shape,'z_scale',z_scale);
    
    % create covariance components that reflect linear shapes
    [~,tmp] = lineImages(100,Dx,k);
    [cc,coeffs] = templates2covariances(tmp,ge.A);
    cc{k+1} = eye(Dx);
    ge.cc = cc;
    
    % create stimuli with omitted gestalt parts
    
end