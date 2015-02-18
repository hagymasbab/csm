function gestalt = gestaltCreate(name,varargin)
    
    % Parse Input
    p = inputParser;
    addParameter(p,'k',2,@isnumeric);
    addParameter(p,'Dx',1024,@isnumeric);
    addParameter(p,'filterShift',2,@isnumeric);
    addParameter(p,'nOrients',4,@isnumeric);
    addParameter(p,'nWLengths',1,@isnumeric);
    addParameter(p,'N',100,@isnumeric);
    addParameter(p,'obsVar',0.1,@isnumeric);
    addParameter(p,'sparsity',0.2,@isnumeric);
    addParameter(p,'z_shape',2,@isnumeric);
    addParameter(p,'z_scale',2,@isnumeric);
    addParameter(p,'g_shape',2,@isnumeric);
    addParameter(p,'g_scale',2,@isnumeric);
    addParameter(p,'null_shape',2,@isnumeric);
    addParameter(p,'null_scale',2,@isnumeric);
    addParameter(p,'B',1,@isnumeric);
    addParameter(p,'filters','line');
    addParameter(p,'prior','gamma');
    addParameter(p,'precision',false,@islogical);
    addParameter(p,'contrast',true,@islogical);
    addParameter(p,'nullComponent',false,@islogical);
    addParameter(p,'overlapping',false,@islogical);
    addParameter(p,'generateComponents',false,@islogical);
    addParameter(p,'componentShape','oriented-gabors');
    addParameter(p,'generateData',false,@islogical);
    addParameter(p,'cc',{});
    parse(p,varargin{:});
    gestalt = p.Results;        
    
    imSize = sqrt(gestalt.Dx); % TODO figure out something if it's not a square
    if strcmp(gestalt.filters,'gabor')        
        %gestalt.A = gaborFilterBank(imSize,imSize,gestalt.filterShift,gestalt.filterShift,[0;pi/4;pi/2;3*pi/4],[4]);
        %gestalt.A = gaborFilterBank(imSize,imSize,gestalt.filterShift,gestalt.filterShift,[pi/4;3*pi/4],[2,4]);
        gestalt.A = gaborFilterBank(imSize,imSize,gestalt.filterShift/2,gestalt.filterShift,[pi/4;3*pi/4],4);
    elseif strcmp(gestalt.filters,'line')        
        gestalt.A = lineFilters(imSize);
    elseif strcmp(gestalt.filters,'eye')
        gestalt.A = eye(gestalt.Dx);
    else
        filterfile = sprintf('filters_%s_%d.mat',gestalt.filters,gestalt.Dx);
        % load filters from file
        load(filterfile);
        gestalt.A = A;
    end
    gestalt.Dv = size(gestalt.A,2);   
    
    % compute some additional matrices to speed up inference
    %fprintf('Calculating A^TA\n');
    %gestalt.A = sparse(gestalt.A);
    gestalt.AA = (gestalt.A)'*gestalt.A;
    %fprintf('Calculating R\n');
    %gestalt.R = pinv(gestalt.AA)*(gestalt.A)';
    %gestalt.tX = zeros(gestalt.N,gestalt.B,gestalt.Dv);
    %for n=1:gestalt.N
    %    gestalt.tX(n,:,:) = reshape(gestalt.X(n,:,:),gestalt.B,gestalt.Dx) * gestalt.A;
    %end
    
    if gestalt.generateComponents
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
            % TODO this is only valid if we have 4-orient 1-frequency
            % gabors with 2-by-2 grid layout
            [gestalt.cc, gestalt.gRF] = gestaltCovariances(gestalt,gestalt.k,'nullComponent',gestalt.nullComponent,'overlapping',gestalt.overlapping,'nOrient',4,'method',gestalt.componentShape);            
    %         if gestalt.nullComponent
    %             gestalt.k = gestalt.k + 1;
    %         end
        end
    elseif isempty(gestalt.cc)
        if gestalt.nullComponent
            gestalt.cc = cell(1,gestalt.k+1);
        else
            gestalt.cc = cell(1,gestalt.k);
        end
    end
    
    gestalt.k = size(gestalt.cc,2);
    gestalt.sparseComponents = length(find(componentSum(1,gestalt.cc))) < gestalt.Dv^2 / 2;
    
    if gestalt.precision
        gestalt.pc = cell(1,gestalt.k);
        for j=1:gestalt.k
            gestalt.pc{j} = inv(gestalt.cc{j});
        end
    end
    
    if gestalt.N > 0 && gestalt.generateData
        gestalt = gestaltGenerate(gestalt,gestalt.N,'batchSize',gestalt.B,'precision',gestalt.precision,'obsVar',gestalt.obsVar,'sparsity',gestalt.sparsity);    
    end
    
    %fprintf('Saving results\n');
    save(strcat('gestalt_',name,'.mat'),'gestalt');