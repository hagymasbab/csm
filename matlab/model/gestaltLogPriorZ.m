function lp = gestaltLogPriorZ(z,ge)
    % this is always a gamma now
    lp = log( gampdf(z,ge.z_shape,ge.z_scale) );
end