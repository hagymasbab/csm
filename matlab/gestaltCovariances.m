function [cc,receptiveFields] = gestaltCovariances(k,Dx,Dv,varargin)
    parser = inputParser;
    addParamValue(parser,'nullComponent',true,@islogical);    
    addParamValue(parser,'overlapping',false,@islogical);    
    parse(parser,varargin{:});        
    params = parser.Results;
    
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
        onindices = zeros(Dv,1);
        if covering
            startpixel = max(g-1,1);
            endpixel = min(g+1,Dv);            
            onindices(startpixel:endpixel) = 1;
        elseif two && ~params.overlapping
            shift = 2;                
            act_shift = 1 + g*shift + (g-1)*width;
            startpixel = (act_shift - 1) * imsizey + vermargin + 1;
            endpixel = act_shift * imsizey - vermargin;
            onindices(startpixel:endpixel) = 1;
        elseif params.overlapping
            if g > 3
                error('Undercomplete representation with more than 3 overlapping components is not implemented.');
            elseif g == 1
                onindices(2*imsizex+2:3*imsizex-2) = 1;
                onindices(3*imsizex+2:4*imsizex-2) = 1;
            elseif g == 2
                for ii = 1:imsizex-4
                    onindices((1+ii)*imsizex+3:(1+ii)*imsizex+4) = 1;
                end
            elseif g == 3
                onindices(3*imsizex+4:4*imsizex-2) = 1;
                onindices(4*imsizex+4:5*imsizex-2) = 1;
            end
        else
            shift = floor((imsizex - k) / (k+1));
            act_shift = 1 + g*shift + (g-1)*width;
            startpixel = (act_shift - 1) * imsizey + vermargin + 1;
            endpixel = act_shift * imsizey - vermargin;
            onindices(startpixel:endpixel) = 1;
        end
        act_rf = zeros(Dv,1);
        act_rf(logical(onindices)) = 1;
        receptiveFields{g} = act_rf;
        
        if params.nullComponent
            cc{g} = 0.01*eye(Dv);
            for pix=1:Dv
                for pix2=pix:Dv
                    if onindices(pix) == 1 && onindices(pix2) == 1
                        cc{g}(pix,pix2) = 1;
                        cc{g}(pix2,pix) = 1;
                    end
                end
                if onindices(pix) == 1
                    cc{g}(pix,pix) = cc{g}(pix,pix) + 1;
                end                
            end
            cc{k+1} = eye(Dv);
        else
            fprintf('....%d/', N);
            idx = find(onindices);
            for nn=1:N
                printCounter(nn);
                x = X(nn,:)';
                x(onindices) = x(idx(1));
                vs(nn,:) = x';
            end
            fprintf('\n....Calculating covariance\n');
            cc{g} = cov(vs);

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
end