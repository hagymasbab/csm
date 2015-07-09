function [max_point,max_ll] = gestaltFindSigmaX(ge,cholesky,X,priorSamples,likeMethod,loadSamples,verbose)
    likefunc = @(sigma) gestaltLogLikelihood2(ge,priorSamples,X,cholesky,'loadSamples',loadSamples,'verbose',0,'method',likeMethod,'sigma',sigma);
    lowerbound = 0.1;
    upperbound = 1;  
    
    res = 10;
    ll = zeros(res,1);
    evalpoints = linspace(lowerbound,upperbound,res);
    for i = 1:res
        if verbose
            printCounter(i,'maxVal',res,'stringVal','Sigma evaluation');
        end
        ll(i) = likefunc(evalpoints(i));
    end
    
    [max_ll,idx] = max(ll);
    max_point = evalpoints(idx);
    
    plot(evalpoints,ll);
    
%     tolerance = 1e-2;
%     smallstep = 1e-2;
%     init_ll = likefunc(sigma_init);
%     max_ll = init_ll;
%     max_point = sigma_init;    
%     
%     % check wheter we should go up or down
%     ll_stepdown = likefunc(sigma_init-smallstep);
%     if ll_stepdown >= init_ll
%         upperbound = sigma_init-smallstep;
%         max_ll = ll_stepdown;
%         max_point = sigma_init-smallstep;
%         upper_ll = ll_stepdown;
%         lower_ll = likefunc(lowerbound);
%     else
%         lowerbound = sigma_init;
%         lower_ll = init_ll;
%         upper_ll = likefunc(upperbound);
%     end
%         
%     while upperbound-lowerbound > tolerance
%         midpoint = (upperbound-lowerbound) / 2;
%         mid_ll = likefunc(midpoint);
%         if mid_
%     end
end