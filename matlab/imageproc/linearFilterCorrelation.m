function cm = linearFilterCorrelation(picNum,patchDBFile,filterBank)
    % load patch database
    load(patchDBFile);
    X = patchDB(:,1:picNum)';
    % load filter bank
    load(filterBank);
    W = inv(A);
    % calculate responses
    U = zeros(picNum,size(X,2));
    for i=1:picNum
        U(i,:) = W * X(i,:)';
    end
    % calculate correlations
    cm = corr(U);
    viewImage(cm);
end
    