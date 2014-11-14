function gestalt = gestaltCreate(name,varargin)
    
    % Parse Input
    p = inputParser;
    addParamValue(p,'k',2,@isnumeric);
    addParamValue(p,'Dx',1024,@isnumeric);
    addParamValue(p,'filterShift',2,@isnumeric);
    addParamValue(p,'nOrients',4,@isnumeric);
    addParamValue(p,'nWLengths',1,@isnumeric);
    addParamValue(p,'N',100,@isnumeric);
    addParamValue(p,'obsVar',0.1,@isnumeric);
    addParamValue(p,'sparsity',0.2,@isnumeric);
    addParamValue(p,'z_shape',1,@isnumeric);
    addParamValue(p,'z_scale',1,@isnumeric);
    addParamValue(p,'g_shape',1,@isnumeric);
    addParamValue(p,'g_scale',0.2,@isnumeric);
    addParamValue(p,'null_shape',2,@isnumeric);
    addParamValue(p,'null_scale',0.2,@isnumeric);
    addParamValue(p,'B',10,@isnumeric);
    addParamValue(p,'filters','line');
    addParamValue(p,'precision',true,@islogical);
    addParamValue(p,'contrast',true,@islogical);
    addParamValue(p,'nullComponent',true,@islogical);
    addParamValue(p,'overlapping',false,@islogical);
    addParamValue(p,'generate',false,@islogical);
    parse(p,varargin{:});
    gestalt = p.Results;        
    
    imSize = sqrt(gestalt.Dx); % TODO figure out something if it's not a square
    if strcmp(gestalt.filters,'gabor')        
        %gestalt.A = gaborFilterBank(imSize,imSize,gestalt.filterShift,gestalt.filterShift,[0;pi/4;pi/2;3*pi/4],[4]);
        gestalt.A = gaborFilterBank(imSize,imSize,gestalt.filterShift,gestalt.filterShift,[pi/4;3*pi/4],[2,4]);
    elseif strcmp(gestalt.filters,'line')        
        gestalt.A = lineFilters(imSize);
    elseif strcmp(gestalt.filters,'eye')
        gestalt.A = eye(gestalt.Dx);
    else
        % load filters from file
        load(gestalt.filters);
        gestalt.A = A;
    end
    gestalt.Dv = size(gestalt.A,2);   
    
    % compute some additional matrices to speed up inference
    %fprintf('Calculating A^TA\n');
    gestalt.AA = (gestalt.A)'*gestalt.A;
    %fprintf('Calculating R\n');
    %gestalt.R = pinv(gestalt.AA)*(gestalt.A)';
    %gestalt.tX = zeros(gestalt.N,gestalt.B,gestalt.Dv);
    %for n=1:gestalt.N
    %    gestalt.tX(n,:,:) = reshape(gestalt.X(n,:,:),gestalt.B,gestalt.Dx) * gestalt.A;
    %end
    
    if gestalt.generate
        if gestalt.Dx == 1
            gestalt.cc{1} = 1;
            gestalt.k = 1;
        elseif gestalt.Dx == 2
            gestalt.cc{1} = [0.9 0.5;0.5 0.3];
            gestalt.k = 1;
        elseif gestalt.Dx == 3
            gestalt.cc{1} = 0.5*[1 0.9 0;0.9 1 0;0 0 1];
            gestalt.cc{2} = 0.5*[1 0 0.9;0 1 0;0.9 0 1];
        elseif gestalt.Dx == 4
            var = 1;
            covar = 1;
            gestalt.cc{1} = var*eye(4) + covar*[0 0 1 0; 0 0 0 0; 1 0 0 0; 0 0 0 0];
            gestalt.cc{2} = var*eye(4) + covar*[0 0 0 0; 0 0 0 1; 0 0 0 0; 0 1 0 0];
        else        
            [gestalt.cc, gestalt.gRF] = gestaltCovariances(gestalt.k,gestalt.Dx,gestalt.Dv,'nullComponent',gestalt.nullComponent,'overlapping',gestalt.overlapping);            
    %         if gestalt.nullComponent
    %             gestalt.k = gestalt.k + 1;
    %         end
        end
    else
        gestalt.cc = cell(1,gestalt.k+1);
    end
    
    gestalt.k = size(gestalt.cc,2);
    
    if gestalt.precision
        gestalt.pc = cell(1,gestalt.k);
        for j=1:gestalt.k
            gestalt.pc{j} = inv(gestalt.cc{j});
        end
    end
    
    if gestalt.generate && gestalt.N > 0
        gestalt = gestaltGenerate(gestalt,gestalt.N,'batchSize',gestalt.B,'precision',gestalt.precision,'obsVar',gestalt.obsVar,'sparsity',gestalt.sparsity);    
    end
    
    %fprintf('Saving results\n');
    save(strcat('gestalt_',name,'.mat'),'gestalt');