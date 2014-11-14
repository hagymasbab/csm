function A = gaborFilterBank(sizex,sizey,dx,dy,orients,wavelengths)
    fprintf('Calculating Gabor filter bank\n');
    if size(orients,1)<size(orients,2)
        orients = orients';
    end
    if size(wavelengths,1)<size(wavelengths,2)
        wavelengths = wavelengths';
    end
    D = sizex * sizey;
    %sigma = floor(sizex/9);
    sigma = 1.5;
    xsteps = floor(sizex/dx);
    ysteps = floor(sizey/dy);
    orNum = size(orients,1);
    freqNum = size(wavelengths,1);
    N = xsteps * ysteps * orNum * freqNum;
    xindexStep = ysteps * orNum * freqNum;
    yindexStep = orNum * freqNum;
    A = zeros(D,N);
    fprintf('..%d/', N);
    for x=1:xsteps
        px = ((x-1)*dx)/sizex;
        for y=1:ysteps
            py = ((y-1)*dy)/sizey;
            for o=1:orNum
                for f=1:freqNum    
                    index = (x-1)*xindexStep + (y-1)*yindexStep + (o-1)*freqNum + f;
                    printCounter(index);
                    %sigma = wavelengths(f) / 2;
                    [~,~,F] = gabor('theta',orients(o),'lambda',wavelengths(f),'width',sizex,'height',sizey,'px',px,'py',py,'Sigma',sigma);         
                    A(:,index) = reshape(F,D,1);
                end
            end
        end
    end
    fprintf('\n');
end