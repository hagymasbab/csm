function gestalt = gestaltCreate(name,varargin)
    
    % Parse Input
    p = inputParser;
    addParamValue(p,'k',2,@isnumeric);
    addParamValue(p,'Dx',32*32,@isnumeric);
    addParamValue(p,'filterShift',2,@isnumeric);
    addParamValue(p,'nOrients',2,@isnumeric);
    addParamValue(p,'nWLengths',2,@isnumeric);
    addParamValue(p,'N',100,@isnumeric);
    addParamValue(p,'obsVar',0.1,@isnumeric);
    addParamValue(p,'sparsity',0.2,@isnumeric);
    addParamValue(p,'B',1,@isnumeric);
    p.KeepUnmatched = true;
    parse(p,varargin{:});
    gestalt = p.Results;        
    
    %gestalt.A = gaborFilterBank(gestalt.imSize,gestalt.imSize,gestalt.filterShift,gestalt.filterShift,[0;pi/2],[4;8]);
    gestalt.A = eye(gestalt.Dx);
    gestalt.Dv = size(gestalt.A,2);   
    
    % compute some additional matrices to speed up inference
    fprintf('Calculating A^TA\n');
    gestalt.AA = (gestalt.A)'*gestalt.A;
    fprintf('Calculating R\n');
    gestalt.R = pinv(gestalt.AA)*(gestalt.A)';
        
    if gestalt.Dx == 1
        gestalt.cc{1} = 1;
        gestalt.k = 1;
    elseif gestalt.Dx == 2
        % TODO
    elseif gestalt.Dx == 3
        gestalt.cc{1} = 0.5*[1 0.9 0;0.9 1 0;0 0 1];
        gestalt.cc{2} = 0.5*[1 0 0.9;0 1 0;0.9 0 1];
    elseif gestalt.Dx == 4
        var = 1;
        covar = 1;
        gestalt.cc{1} = var*eye(4) + covar*[0 0 1 0; 0 0 0 0; 1 0 0 0; 0 0 0 0];
        gestalt.cc{2} = var*eye(4) + covar*[0 0 0 0; 0 0 0 1; 0 0 0 0; 0 1 0 0];
    else        
        gestalt.cc = gestaltCovariances(gestalt.k,gestalt.R);
    end
    
    gestalt = gestaltGenerate(gestalt,gestalt.N);
    fprintf('Transforming synthetic data\n');
    gestalt.tX = zeros(gestalt.N,gestalt.B,gestalt.Dv);
    for n=1:gestalt.N
        gestalt.tX(n,:,:) = reshape(gestalt.X(n,:,:),gestalt.B,gestalt.Dx) * gestalt.A';
    end
    
    fprintf('Saving results\n');
    save(strcat('gestalt_',name,'.mat'),'gestalt');