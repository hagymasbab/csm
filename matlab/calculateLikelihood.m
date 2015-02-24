function calculateLikelihood(iterfile,gefile,dataset,N_test,samples,randseed,calculate)

    setrandseed(randseed);
    namepart = iterfile(5:end);
    
    if strcmp(calculate,'likelihood')
        load(dataset);
        N_all = size(patchDB,2);
        X_test = patchDB(:,chooseKfromN(N_test,N_all))';        
    end
    
    load(iterfile);
    load(gefile);
    ge = ge_saved;
    
    % find out how many iterations are actually there in the data
    calc = true;
    index = 1;
    while calc
        if isempty(state_sequence{index})
            calc = false;
            index = index - 1;
            continue;
        else            
            index = index + 1;
        end
    end
    
    
    likelihoods = zeros(index,1);
    norms = zeros(index,ge.k);
    distances = zeros(index,ge.k*(ge.k-1)/2);
    
    for i=1:index
        printCounter(i,'stringVal','Iteration','maxVal',index);
        act_cc = state_sequence{i}.estimated_components;
        if strcmp(calculate,'likelihood')
            ge.cc = act_cc;
            ll = gestaltLogLikelihood(ge,samples,X_test,'scientific',true);        
            likelihoods(i) = ll;
            %save('iter_likes.mat','likelihoods');
            save(['bin/likes' namepart],'likelihoods');
        elseif strcmp(calculate,'metrics')
            met_idx = 1;
            for kk = 1:ge.k
              norms(i,kk) = norm(act_cc{kk});
              for other = kk+1:ge.k
                  dist = covcompRootMeanSquare(act_cc(kk),act_cc(other),1);
                  distances(i,met_idx) = dist;
                  met_idx = met_idx+1;
              end
            end
            %save('iter_norms.mat','norms','distances');
            save(['bin/norms' namepart],'norms','distances');
        elseif strcmp(calculate,'grf')            
            if i == 1
                pixel_comps = {};
                angular_stds = zeros(index,ge.k,3);
                pixel_rms = zeros(index-1,ge.k);
            end
            [~,seeds,angstds,~] = gestaltGReceptiveFields(ge,act_cc,10000,false);
            pixel_comps{end+1} = seeds;
            angular_stds(i,:,:) = angstds;
            if i>1
                for j = 1:ge.k
                    pixel_rms(i-1,j) = sqrt(sum((pixel_comps{i-1}{j} - pixel_comps{i}{j}).^2));
                end
            end
            %save('iter_grf.mat','pixel_comps','angular_stds','pixel_rms');
            save(['bin/grf' namepart],'pixel_comps','angular_stds','pixel_rms');
        end
    end
end