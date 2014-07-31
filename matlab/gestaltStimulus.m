function X = gestaltStimulus(Dx,B,left,half)
    imd = sqrt(Dx);
    X = zeros(B,Dx);
    for b=1:B
        act = rand(imd);
        horcoord = 3;
        if ~left
            horcoord = imd-2;
        end
        verend = imd - 1;
        if half
            verend = floor(imd/2);
        end
        act(2:verend,horcoord) = 2;
        X(b,:) = act(:)';
    end
end