function [cc,receptiveFields] = gestaltCovariances(ge,k,varargin)
    parser = inputParser;
    addParameter(parser,'nullComponent',true,@islogical);    
    addParameter(parser,'overlapping',false,@islogical);    
    addParameter(parser,'verbose',false,@islogical);    
    addParameter(parser,'nOrient',4,@isnumeric);    
    addParameter(parser,'method','oriented-gabors');  
    parse(parser,varargin{:});        
    params = parser.Results;
    
    covering = false;
    two = false;
    if k == 0;
        % we will generate a component set that covers the whole field in a
        % way that every V-unit belongs to 2 vertical components of length
        % 4 or 2 on the edges
        k = ge.Dv;
        covering = true;
    elseif k == 2
        two = true;
    end                
    
    if params.verbose
        fprintf('Calculating covariance components\n');
    end   
   
    if strcmp(params.method,'oriented-gabors')
        cc = {};
        sh1 = 1; % single shift 
        sh2 = 2; % double shift
        imdim = sqrt(ge.Dx);                  
        shift = params.nOrient / 2;
        nRF = ge.Dx / params.nOrient;
        rfsinarow = imdim/shift;        
        
        rfs = 1:nRF;
        last_orients = 4;
        if ~covering
            used_rfs = ceil(k/params.nOrient);
            rfs = chooseKfromN(used_rfs,nRF);
            last_orients = k - (params.nOrient*(used_rfs-1));
        end
        
        for i = 1:length(rfs);
            ic_idx = rfs(i);
            % define gestalts as sets of RF locations and orientation indices
            % this only works with the standard 4-oriented gabors
            templates = [ic_idx-sh2                 ic_idx-sh1                 ic_idx+sh1                 ic_idx+sh2;                 ...
                         ic_idx-(sh2*rfsinarow)+sh2 ic_idx-(sh1*rfsinarow)+sh1 ic_idx+(sh1*rfsinarow)-sh1 ic_idx+(sh2*rfsinarow)-sh2; ...
                         ic_idx-(sh2*rfsinarow)     ic_idx-(sh1*rfsinarow)     ic_idx+(sh1*rfsinarow)     ic_idx+(sh2*rfsinarow);     ...
                         ic_idx-(sh2*rfsinarow)-sh2 ic_idx-(sh1*rfsinarow)-sh1 ic_idx+(sh1*rfsinarow)+sh1 ic_idx+(sh2*rfsinarow)+sh2];
            if i == length(rfs)
                templates = templates(1:last_orients,:);
            end
            templates = max(templates,1);
            templates = min(templates,nRF);
            for o=1:size(templates,1)
                gestaltvec = zeros(ge.Dv,1);
                filter_indices = (templates(o,:)' - 1) * params.nOrient + o;
                gestaltvec(filter_indices) = 1;
                act_cc = gestaltvec*gestaltvec' * 10;
                if ~params.nullComponent
                    %act_cc = act_cc + (1/ge.Dv) * eye(ge.Dv);
                    act_cc = act_cc + (1/sqrt(ge.Dv)) * eye(ge.Dv);
                end
                cc{end+1} = act_cc;
            end
        end
        if params.nullComponent
            cc{end+1} = eye(ge.Dv);
        end
        ge.k = length(cc);
        % calculate the receptive fields too
        % receptiveFields = gestaltGReceptiveFields(ge,cc,1000,false);        
    
    elseif strcmp(params.method,'block')
        if params.overlapping            
            blocksize = floor(1.5*(ge.Dx-k+1) / k);
            starts = ((1:k)-1)*(floor(blocksize/2)+1) + 2;            
        else
            blocksize = floor((ge.Dx-k+1) / k);
            starts = ((1:k)-1)*(blocksize+1) + 2;            
        end
        
        ends = starts + blocksize -1;            
        cc = {};
        for i=1:k
            act_cc = zeros(ge.Dv);
            act_cc(starts(i):ends(i),starts(i):ends(i)) = 2;
            diagonals = sum(act_cc);
            act_cc = act_cc + diag(diagonals);
            cc{i} = act_cc;
        end
        
        if params.nullComponent
            cc{end+1} = eye(ge.Dv);
        end
        ge.k = length(cc);
        
        receptiveFields = gestaltGReceptiveFields(ge,cc,1000,false);   
        
    elseif strcmp(params.method,'vertical-bars')
        
         % vertical lines
        imsizex = floor(sqrt(ge.Dx));
        imsizey = ceil(sqrt(ge.Dx));

        width = 1;    
        vermargin = 1;

        N = max(ge.Dx,ge.Dv) + 1;
        for g = 1:k
            if params.verbose
                fprintf('..Component %d\n', g);
            end

            vs = zeros(N,ge.Dv);
            X = mvnrnd(zeros(N,ge.Dx),eye(ge.Dx));
            onindices = zeros(ge.Dv,1);
            if covering
                startpixel = max(g-1,1);
                endpixel = min(g+1,ge.Dv);            
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
            act_rf = zeros(ge.Dv,1);
            act_rf(logical(onindices)) = 1;
            receptiveFields{g} = act_rf;

            if params.nullComponent
                cc{g} = 0.01*eye(ge.Dv);
                for pix=1:ge.Dv
                    for pix2=pix:ge.Dv
                        if onindices(pix) == 1 && onindices(pix2) == 1
                            cc{g}(pix,pix2) = 1;
                            cc{g}(pix2,pix) = 1;
                        end
                    end
                    if onindices(pix) == 1
                        cc{g}(pix,pix) = cc{g}(pix,pix) + 1;
                    end                
                end
                cc{k+1} = eye(ge.Dv);
            else
                fprintf('....%d/', N);
                igeDx = find(onindices);
                for nn=1:N
                    printCounter(nn);
                    x = X(nn,:)';                    
                    x(logical(onindices)) = x(igeDx(1));
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
                    cc{g} = cc{g} + eyecoeff * eye(ge.Dv);
                end
            end
        end
    else
        error('nonexistent cc method');
    end
end