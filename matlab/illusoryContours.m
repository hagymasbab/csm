function illusoryContours(randseed,nTrials,nSamples,nCont,pre_cs,nullComp,backZ)
    %close all;
    Dx = 1024;    
    %nullComp = false;

    if strcmp(randseed,'last')
        load lastrandseed;
    end
    s = RandStream('mt19937ar','Seed',randseed);
    RandStream.setGlobalStream(s);
    randseed = s.Seed;
    save('lastrandseed.mat','randseed');
    
    [stimuli,A,gestalts] = ICstimuli(nCont,2,backZ);
    
    % using fewer stimuli
    numstim = 3;
    stimuli = stimuli(2:3,1:numstim);
    gestalts = gestalts(1:numstim,:);
    
    save('ic.mat','A');
    k = size(gestalts,1);
    central_field = gestalts(:,1);
    activated_units = gestalts(:,2:end);
    activated_units = activated_units(:)';

%     %A = gaborFilterBank(32,32,1,2,[pi/4,3*pi/4],[4]);
%     central_cell = 8*64+17;
%     upper_cell = 12*64+9;
%     omitted_cell = 10*64+13;
%     lower_cell = 4*64+25;
%     other_cell = 6*64+21;
%     gestalts = [central_cell upper_cell omitted_cell;central_cell lower_cell other_cell];
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
        if ~nullComp
            cc{kk} = cc{kk} + 0.1*eye(Dx);
        end
    end
    if nullComp
        cc{end+1} = eye(Dx);
    end
    
    % model parameters for generation and sampling    
    generating_sigma = 0.1;
    z_gen = 1;
    g_gen = 10;
        
    sampling_sigma = 1;
    g_shape = 1;
    g_scale = 2;
    z_shape = 1;
    z_scale = 0.1;
    sample_z = true;
    initZ = 10; %in case we don't sample, this remains the same
    filter_scaling = 1;
    g_init = zeros(k,1);
    if nullComp
        g_init = [g_init;1];
    else
        g_init = g_init + 0.1;
        %g_init = [];
    end
    
    ge = gestaltCreate('temp','Dx',Dx,'k',k,'B',1,'N',1,'nullComponent',nullComp, ...
        'filters','ic.mat','obsVar',generating_sigma,'g_shape',g_shape,'g_scale',g_scale,'z_shape',z_shape,'z_scale',z_scale);
    ge.cc = cc;
    
    prestimSamp = 2;
    poststimSamp = 3;
%     g_stim = [g_gen;0;0];
%     [ic_stim,gen_v] = gestaltAncestralSample(ge,g_stim,z_gen,false,false);
%     ic_stim = ic_stim' - gen_v(1,omitted_cell) * ge.A(:,omitted_cell);
%     viewImage(ic_stim);
    
%     A = ge.A;
%     save('ic.mat','A','ic_stim','gestalts');
    
    nStim = size(stimuli,2);
    nCont = size(stimuli,1);
    %nCont = size(stimuli,1);
    ge.obsVar = sampling_sigma;
    
    allsamp = zeros(nStim,nCont,nTrials,nSamples,ge.k+ge.Dv);
    zsamp = zeros(nStim,nCont,nTrials,nSamples);
    within_trial_variance = zeros(nStim,nCont,nTrials,ge.k+ge.Dv);
    trial_to_trial_variance = zeros(nStim,nCont,nSamples,ge.k+ge.Dv);
    within_trial_covariance = zeros(nStim,nCont,nTrials,ge.k+ge.Dv,ge.k+ge.Dv);
    for stim = 1:nStim
        for cont = 1:nCont
            fprintf('Stimulus %d/%d, contrast %d/%d ',nStim,stim,nCont,cont);
            ge.X(1,:,:) = reshape(stimuli{cont,stim},1,Dx);
            
%             central_samples = zeros(nTrials,nSamples);
%             omitted_samples = zeros(nTrials,nSamples);
%             other_samples = zeros(nTrials,nSamples);
%             g1_samples = zeros(nTrials,nSamples);
%             g2_samples = zeros(nTrials,nSamples);
%             final_v = zeros(nTrials,ge.Dv);
            
            for t = 1:nTrials
                printCounter(t,'stringVal','trial','maxVal',nTrials,'newLine',true);
                if isempty(pre_cs)
                    [cs,~,zs] = gestaltGibbs(ge,1,nSamples,'verbose',0,'vSampler','direct','contrast',sample_z, ...
                            'prestimSamples',prestimSamp,'poststimSamples',poststimSamp,'verbose',1,'initZ',initZ,'initG',g_init);
                else
                    cs = squeeze(pre_cs(stim,cont,t,:,:));
                end

%                 central_samples(t,:) = cs(:,ge.k + central_cell)';
%                 omitted_samples(t,:) = cs(:,ge.k + omitted_cell)';
%                 other_samples(t,:) = cs(:,ge.k + other_cell)';
%                 g1_samples(t,:) = cs(:,1)';
%                 g2_samples(t,:) = cs(:,2)';
%                 final_v(t,:) = cs(end-poststimSamp-1,ge.k+1:end);
                allsamp(stim,cont,t,:,:) = cs;
                zsamp(stim,cont,t,:) = zs;
                within_trial_variance(stim,cont,t,:) = var(cs,0,1);
                within_trial_covariance(stim,cont,t,:,:) = cov(cs);                
%                 fprintf('\n');
            end
            trial_to_trial_variance(stim,cont,:,:) = var(allsamp(stim,cont,:,:,:),0,3);
            save('ic_samp.mat','allsamp','zsamp','trial_to_trial_variance','within_trial_variance','within_trial_covariance', ...
                'central_field','A','k','prestimSamp','poststimSamp','activated_units','nullComp');
            
%             figure();
%             row = 3;
%             col = 2;
%             plotPair(central_samples,omitted_samples,row,col,1,true,{'stimulated v','gestalt activated v'},{prestimSamp,poststimSamp},{});
%             plotPair(other_samples,omitted_samples,row,col,3,true,{'other gestalt v','gestalt activated v'},{prestimSamp,poststimSamp},{});
%             plotPair(g1_samples,g2_samples,row,col,5,true,{'activated g','non-activated g'},{prestimSamp,poststimSamp},{});
% 
%             figure();
%             final_percept = ge.A * mean(final_v)';
%             viewImage(final_percept);
        end
    end
    
end
    