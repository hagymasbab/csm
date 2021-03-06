function X = gestaltStimulus(Dx,B,left,half,addIrrel)
    imd = sqrt(Dx);
    X = zeros(B,Dx);
    randvals = randn(B,Dx);
    for b=1:B
        act = reshape(randvals(b,:),imd,imd);
        horcoord = 3;        
        verend = imd - 1;
        if half
            verend = floor(imd-3);
        end
        act(2:verend,horcoord) = 2 + randn(1);
        %act(2:verend,horcoord) = act(2:verend,horcoord) + 1;
        %act(2:verend,horcoord) = act(2,horcoord);
        
        if addIrrel
            act(2:verend,horcoord+1) = 2 + randn(1);
        end
        
        if ~left
            horcoord = imd-2;
            act(2:verend,horcoord) = 2 + randn(1);
        end
        
        X(b,:) = act(:)';
    end
end