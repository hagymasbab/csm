function [gsamp,zsamp] = gestaltRejectionGZ(x,ge,L,randseed)
    % find the MAP
    options = optimoptions('fmincon','Algorithm','interior-point','Display', 'off');
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
    vars_opt = fmincon(targetfunc,init,constraints_A,constraints_b,[],[],[],[],[],options)
end