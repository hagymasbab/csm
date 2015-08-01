function [cm,U] = linearFilterCorrelation(picNum,patchDBFile,filterBank,gsm_res)
    % load patch database
    load(patchDBFile);
    X = patchDB(:,1:picNum)';
    % load filter bank
    load(filterBank);
    %W = pinv(A);
    W = A';
    % calculate responses
    U = zeros(picNum,size(X,2));
    for i=1:picNum
        printCounter(i,'maxVal',picNum,'stringVal','Patch');
        if gsm_res == 0
            U(i,:) = W * X(i,:)';
        else
            U(i,:) = gsmPosteriorV(X(i,:)',A,eye(size(A,2)),0.5,2,2,gsm_res)';
        end
    end
    % calculate correlations
    cm = corr(U);
    %viewImage(cm);
    %cm(856,980)
    %scatter(U(:,856),U(:,980));
end
    