function img = grating(lambda,theta,phase,imSize)
    % imSize/2 >= lambda >= 3
    % 180 > theta >= 0
    % 1 > phase >= 0
    freq = imSize/lambda;
    phaseRad = (phase * 2* pi);
    X0 = 1:imSize;
    X0 = (X0 / imSize) - .5;
    [Xm,Ym] = meshgrid(X0, X0); 
    thetaRad = (theta / 360) * 2*pi;        % convert theta (orientation) to radians
    Xt = Xm * cos(thetaRad);                % compute proportion of Xm for given orientation
    Yt = Ym * sin(thetaRad);                % compute proportion of Ym for given orientation
    XYt = Xt + Yt;                          % sum X and Y components
    XYf = XYt * freq * 2*pi;                % convert to radians and scale by frequency
    img = sin( XYf + phaseRad);    
end