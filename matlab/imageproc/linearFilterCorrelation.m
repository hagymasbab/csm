function [cm,U] = linearFilterCorrelation(picNum,patchDBFile,filterBank)
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
        U(i,:) = W * X(i,:)';
    end
    % calculate correlations
    cm = corr(U);
    %viewImage(cm);
    %cm(856,980)
    %scatter(U(:,856),U(:,980));
end
    