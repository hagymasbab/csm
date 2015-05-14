function cholesky = cellchol(cc,varargin)
    parser = inputParser;   
    addParameter(parser,'cheat',false,@islogical); 
    parse(parser,varargin{:});        
    params = parser.Results;  
    
    cholesky = cc;
    for j=1:size(cc,2)
        if params.cheat 
            [~,e] = chol(cholesky{j});
            if e~= 0
                cholesky{j} = nearestSPD(cholesky{j});
            end
        end
        cholesky{j} = chol(cholesky{j});
    end    
end