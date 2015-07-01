function gestaltFindSigmaX(ge,cholesky,sigma_init,X)
    likefunc = @(sigma) gestaltLogLikelihood2(ge,params.priorSamples,X,cholesky,'loadSamples',true,'verbose',0,'method','intuition','sigma',sigma);
    % check wheter we should go up or down
    % check some point far away, if it's still better, find a point where
    % it's not
    % recursively look for the extreme value
end