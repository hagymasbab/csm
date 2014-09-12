function [cc,receptiveFields] = gestaltCovariances(k,Dx,Dv,nullComponent)
    covering = false;
    two = false;
    if k == 0;
        % we will generate a component set that covers the whole field in a
        % way that every V-unit belongs to 2 vertical components of length
        % 4 or 2 on the edges
        k = Dv;
        covering = true;
    elseif k == 2
        two = true;
    end                

    fprintf('Calculating covariance components\n');
    % vertical lines
    imsizex = floor(sqrt(Dx));
    imsizey = ceil(sqrt(Dx));
    
    width = 1;    
    vermargin = 1;
   
    N = max(Dx,Dv) + 1;
    for g = 1:k
        fprintf('..Component %d\n', g);
                
        vs = zeros(N,Dv);
        X = mvnrnd(zeros(N,Dx),eye(Dx));
        if covering
            startpixel = max(g-1,1);
            endpixel = min(g+1,Dv);            
        elseif two
            shift = 2;                
            act_shift = 1 + g*shift + (g-1)*width;
            startpixel = (act_shift - 1) * imsizey + vermargin + 1;
            endpixel = act_shift * imsizey - vermargin;
        else
            shift = floor((imsizex - k) / (k+1));
            act_shift = 1 + g*shift + (g-1)*width;
            startpixel = (act_shift - 1) * imsizey + vermargin + 1;
            endpixel = act_shift * imsizey - vermargin;
        end
        act_rf = zeros(Dv,1);
        act_rf(startpixel:endpixel) = 1;
        receptiveFields{g} = act_rf;
        
        if nullComponent
            cc{g} = 0.01*eye(Dv);
            length = endpixel-startpixel + 1;
            cc{g}(startpixel:endpixel,startpixel:endpixel) = ones(length) + eye(length);
            continue;
        end
        fprintf('....%d/', N);
        for i=1:N
            printCounter(i);
            x = X(i,:);
            x(1,startpixel:endpixel) = x(1,startpixel);
%             if covering                
%                 x(1,startpixel:endpixel) = x(1,startpixel);
%             else
%                 x = reshape(x,imsizex,imsizey);
%                 x(vermargin+1:imsizex-vermargin,act_shift:act_shift+width-1) = x(vermargin+1,act_shift);
%                 x = reshape(x,1,Dx);
%             end
                                    
            %vs(i,:) = x * R';      
            vs(i,:) = x;
        end
        fprintf('\n....Calculating covariance\n');
        cc{g} = cov(vs);
        
        if nullComponent
            for pix=1:Dv
                if pix < startpixel || pix > endpixel
                    cc{g}(pix,pix) = 0.01;
                end
            end
        else
            % preventing ill-conditioned covariance components
            if rcond(cc{g}) < 1e-10
                fprintf('ill-conditioned\n');
                maxelem = max(cc{g}(:));
                mindiag = min(diag(cc{g}));
                %eyecoeff = maxelem-mindiag;
                eyecoeff = 0.1;
                cc{g} = cc{g} + eyecoeff * eye(Dv);
            end
        end
    end
    
    if nullComponent
        cc{k+1} = eye(Dv);
    end
end