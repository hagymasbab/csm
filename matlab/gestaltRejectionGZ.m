function [gsamp,zsamp] = gestaltRejectionGZ(x,ge,L,randseed)
    setrandseed(randseed);

    fprintf('Defining the proposal...');
    % find the MAP
    options = optimoptions('fmincon','Algorithm','sqp','Display', 'off');
    %options = optimoptions('fmincon','Algorithm','sqp');
    constraints_A = [];
    constraints_b = [];    
    % all g and z should be positive or zero
    for kk = 1:ge.k+1
        act_A = zeros(1,ge.k+1);
        act_A(kk) = -1;
        constraints_A = [constraints_A; act_A];
        constraints_b = [constraints_b; 0];
    end
    % the sum of g should be nonzero
    constraints_A = [constraints_A; [-ones(1,ge.k) 0]];
    constraints_b = [constraints_b; 0.1];
    
    init = [ge.g_scale*ones(ge.k,1);ge.z_scale];
    targetfunc = @(vars) -gestaltLogPostGZ(vars(1:ge.k,1),vars(end,1),x,ge); 
    [vars_opt, target_at_map] = fmincon(targetfunc,init,constraints_A,constraints_b,[],[],[],[],[],options);
        
%     vars_opt
    target_at_map = -target_at_map
    
    % define the constant for a Gamma proposal
    alpha = 2;
    theta = vars_opt;
    proposalfunc = @(vars) sum(log(gampdf(vars,alpha,theta)));
    
    proposal_at_map = proposalfunc(vars_opt)
    log_konst = target_at_map - proposal_at_map;
    
    target_at_10 = -targetfunc(vars_opt+10)
    proposal_at_10 = proposalfunc(vars_opt+10)
    difference_at_10 = target_at_10 - proposal_at_10;
    target_at_001 = -targetfunc(0.00001*ones(ge.k+1,1))
    proposal_at_001 = proposalfunc(0.00001*ones(ge.k+1,1))
    difference_at_001 = target_at_001 - proposal_at_001;
    
%     fprintf('done. ');
% 
%     gsamp = zeros(L,ge.k);
%     zsamp = zeros(L,1);
%     rejections = 0;
%     count = 0;
%     while count < L        
%         proposition = gamrnd(alpha,theta);
%         target_at_proposition = -targetfunc(proposition)
%         proposal_at_proposition = log_konst + proposalfunc(proposition)
%         upperbound = exp(proposal_at_proposition)
%         auxiliary = upperbound * rand(1)
%         log(auxiliary)
%         pause
%         if log(auxiliary) <= target_at_proposition
%             count = count + 1;
%             printCounter(count,'maxVal',L,'stringVal','Sample');
%             gsamp(count,:) = proposition(1:ge.k,1);
%             zsamp(count,:) = proposition(end,1);
%         else
%             rejections = rejections + 1;
%         end
%     end
%             
%     rejections
%     mean(gsamp)
%     mean(zsamp)
    
end