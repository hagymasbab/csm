function cc = gestaltCovariances(k,R,plusShift)
    fprintf('Calculating covariance components\n');
    % vertical lines
    Dx = size(R,2);
    imsizex = floor(sqrt(Dx));
    imsizey = ceil(sqrt(Dx));
    Dv = size(R,1);
    %fprintf('Calculating R\n');
    % R = pinv(A'*A)*A';
    % the gestalts should be placed over or watermarked onto natural images
    shift = floor(imsizex/(k+1));
    width = max(1,floor(shift/2));
    margin = width;
    N = max(Dx,Dv) + 1;
    for g = 1:k
        fprintf('..Component %d\n', g);
        act_shift = g*shift + plusShift;
        vs = zeros(N,Dv);
        X = mvnrnd(zeros(N,Dx),eye(Dx));
        fprintf('....%d/', N);
        for i=1:N
            printCounter(i);
            x = reshape(X(i,:),imsizex,imsizey);
            x(margin:imsizex-margin,act_shift:act_shift+width) = x(margin,act_shift);
            x = reshape(x,1,Dx);
            %vs(i,:) = x * R';      
            vs(i,:) = x;
        end
        fprintf('\n....Calculating covariance\n');
        cc{g} = cov(vs);
        % preventing ill-conditioned covariance components
        if rcond(cc{g}) < 1e-10
            cc{g} = cc{g} + eye(Dv);
        end
    end
end