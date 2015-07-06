function [expectation,weights] = importanceSampling(func,target_pdf,proposal_pdf,proposal_samp,L)
    fprintf('untested');
    
    % take L samples from the proposal
    prop_samp = proposal_samp(L);
    
    % calculate the importance weights from target values
    func_eval = zeros(L,1);
    weights = zeros(L,1);
    summed_weights = 0;
    for i=1:L        
        act_samp = prop_samp(i,:)';
        act_target = target_pdf(act_samp);
        act_prop = proposal_pdf(act_samp);
        weights(i) = act_target / act_prop;
        summed_weights = summed_weights + weights(i);
        func_eval(i) = func(act_samp);
    end
    weights = weights / summed_weights;
    
    % calculate the expectation
    expectation = weights' * func_eval;
    
end