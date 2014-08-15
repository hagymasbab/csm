function cc = gestaltCovariances(k,Dx,Dv)
    covering = false;
    if k == Dv/2 + 1;
        % we will generate a component set that covers the whole field in a
        % way that every V-unit belongs to 2 vertical components of length
        % 4 or 2 on the edges
        covering = true;
    end                

    fprintf('Calculating covariance components\n');
    % vertical lines
    imsizex = floor(sqrt(Dx));
    imsizey = ceil(sqrt(Dx));
    %fprintf('Calculating R\n');
    % R = pinv(A'*A)*A';
    % the gestalts should be placed over or watermarked onto natural images
    %shift = floor(imsizex/(k+1));
    shift = 2;
    %width = max(1,floor(shift/2));
    width = 1;
    %margin = width;
    vermargin = 1;
    N = max(Dx,Dv) + 1;
    for g = 1:k
        fprintf('..Component %d\n', g);
        act_shift = 1 + g*shift + (g-1)*width;
        vs = zeros(N,Dv);
        X = mvnrnd(zeros(N,Dx),eye(Dx));
        fprintf('....%d/', N);
        for i=1:N
            printCounter(i);
            x = X(i,:);
            if covering
                startpixel = max((g-2) * 2,1);
                endpixel = min((g-2) * 2+4,Dx);
                x(1,startpixel:endpixel) = x(1,startpixel);
            else
                x = reshape(x,imsizex,imsizey);
                x(vermargin+1:imsizex-vermargin,act_shift:act_shift+width-1) = x(vermargin+1,act_shift);
                x = reshape(x,1,Dx);
            end
                                    
            %vs(i,:) = x * R';      
            vs(i,:) = x;
        end
        fprintf('\n....Calculating covariance\n');
        cc{g} = cov(vs);
        % preventing ill-conditioned covariance components
        if rcond(cc{g}) < 1e-10
            fprintf('\nill-conditioned\n');
            maxelem = max(cc{g}(:));
            mindiag = min(diag(cc{g}));
            %eyecoeff = maxelem-mindiag;
            eyecoeff = 0.1;
            cc{g} = cc{g} + eyecoeff * eye(Dv);
        end
    end
end