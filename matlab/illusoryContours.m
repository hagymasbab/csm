function illusoryContours(randseed,nTrials)
    close all;
    Dx = 1024;
    k = 2;

    if strcmp(randseed,'last')
        load lastrandseed;
    end
    s = RandStream('mt19937ar','Seed',randseed);
    RandStream.setGlobalStream(s);
    randseed = s.Seed;
    save('lastrandseed.mat','randseed');

    %A = gaborFilterBank(32,32,1,2,[pi/4,3*pi/4],[4]);
    central_cell = 8*64+17;
    upper_cell = 12*64+9;
    omitted_cell = 10*64+13;
    lower_cell = 4*64+25;
    other_cell = 6*64+21;
    gestalts = [central_cell upper_cell omitted_cell;central_cell lower_cell other_cell];
    cc = {};
    for kk=1:k
        cc{kk} = zeros(Dx);
        for i = 1:3
            cc{kk}(gestalts(kk,i),gestalts(kk,i)) = 1;
            for j =i+1:3
               cc{kk}(gestalts(kk,i),gestalts(kk,j)) = 0.5; 
               cc{kk}(gestalts(kk,j),gestalts(kk,i)) = cc{kk}(gestalts(kk,i),gestalts(kk,j));
            end
        end
    end
    cc{end+1} = eye(Dx);
    
    % model parameters for generation and sampling
    generating_sigma = 0.1;
    z_gen = 1;
    g_gen = 10;
        
    sampling_sigma = 1;
    g_scale = 2;
    z_shape = 1;
    z_scale = 0.1;
    sample_z = true;
    fixedZ = 0.1; %in case we don't sample  
    filter_scaling = 1;
    g_init = [zeros(k,1);1];
    
    ge = gestaltCreate('temp','Dx',Dx,'k',k,'B',1,'N',1, ...
        'filters','gabor_1024.mat','obsVar',generating_sigma,'g_scale',g_scale,'z_shape',z_shape,'z_scale',z_scale);
    ge.cc = cc;
    
    nSamples = 13;
    prestimSamp = 3;
    g_stim = [g_gen;0;0];
    [ic_stim,gen_v] = gestaltAncestralSample(ge,g_stim,z_gen,false,false);
    ic_stim = ic_stim' - gen_v(1,omitted_cell) * ge.A(:,omitted_cell);
    viewImage(ic_stim);
    
    A = ge.A;
    save('ic.mat','A','ic_stim','gestalts');
    
    ge.X(1,:,:) = ic_stim';
    ge.obsVar = sampling_sigma;
    central_samples = zeros(nTrials,nSamples);
    omitted_samples = zeros(nTrials,nSamples);
    other_samples = zeros(nTrials,nSamples);
    g1_samples = zeros(nTrials,nSamples);
    g2_samples = zeros(nTrials,nSamples);
    final_v = zeros(nTrials,ge.Dv);
    allsamp = zeros(nTrials,nSamples,ge.k+ge.Dv);
    for t = 1:nTrials
        [cs,~,cz] = gestaltGibbs(ge,1,nSamples,'verbose',0,'vSampler','direct','contrast',sample_z, ...
                    'prestimSamples',prestimSamp,'verbose',1,'fixedZ',fixedZ,'initG',g_init);
        central_samples(t,:) = cs(:,ge.k + central_cell)';
        omitted_samples(t,:) = cs(:,ge.k + omitted_cell)';
        other_samples(t,:) = cs(:,ge.k + other_cell)';
        g1_samples(t,:) = cs(:,1)';
        g2_samples(t,:) = cs(:,2)';
        final_v(t,:) = cs(end,ge.k+1:end);
        allsamp(t,:,:) = cs;
        save('ic_samp.mat',allsamp);
        fprintf('\n');
    end
    figure();
    row = 3;
    col = 2;
    plotPair(central_samples,omitted_samples,row,col,1,true,{'stimulated v','gestalt activated v'},prestimSamp,{});
    plotPair(other_samples,omitted_samples,row,col,3,true,{'other gestalt v','gestalt activated v'},prestimSamp,{});
    plotPair(g1_samples,g2_samples,row,col,5,true,{'activated g','non-activated g'},prestimSamp,{});
    
    figure();
    final_percept = ge.A * mean(final_v)';
    viewImage(final_percept);
    
end
    