function calculateLikelihood(iterfile,gefile,dataset,N_test,samples,randseed,calculate)

    setrandseed(randseed);
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
            save('iter_likes.mat','likelihoods');
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
            save('iter_norms.mat','norms','distances');
        elseif strcmp(calculate,'grf')
            if i>1
%                 [~,seeds,angstds,~] = gestaltGReceptiveFields(ge,act_cc,10000,false);
            else
%                 angstds = circ_rad2ang(sqrt(2)) * ones(ge.k,3);
%                 seeds = {};
%                 for kk = 1:ge.k
%                     seeds{end+1} = zeros(ge.Dv,1);
%                 end
                pixel_comps = {};
                angular_stds = zeros(index,ge.k,3);
            end
            [~,seeds,angstds,~] = gestaltGReceptiveFields(ge,act_cc,10000,false);
            pixel_comps{end+1} = seeds;
            angular_stds(i,:,:) = angstds;
            save('iter_grf.mat','pixel_comps','angular_stds');
        end
    end
end