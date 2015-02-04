function calculateLikelihood(iterfile,gefile,dataset,N_test,samples,randseed)

    setrandseed(randseed);

    load(dataset);
    N_all = size(patchDB,2);
    X_test = patchDB(:,chooseKfromN(N_test,N_all))';        
    
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
    
    for i=1:index
        printCounter(i,'stringVal','Iteration','maxVal',46);
        act_cc = state_sequence{i}.estimated_components;
        ge.cc = act_cc;
        ll = gestaltLogLikelihood(ge,samples,X_test,'scientific',true);        
        likelihoods(i) = ll;
        save('iter_likes.mat','likelihoods');
    end
end