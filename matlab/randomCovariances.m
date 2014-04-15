function cc = randomCovariances(k,dim,varargin)
    p = inputParser;
    addParamValue(p,'blockStrength',0,@isnumeric);
    addParamValue(p,'precision',false,@islogical);
    parse(p,varargin{:});
    blockStrength = p.Results.blockStrength;
    precision = p.Results.precision;

    blocksize = floor(dim/k);
    block = blockStrength * ones(blocksize);
    for i=1:k
          base = rand(dim,dim);
          cc{i} = (0.1/dim)*(base'*base);
          cc{i} = cc{i} + eye(dim);
          begin_block = (i-1)*blocksize + 1;
          end_block = i*blocksize;          
          cc{i}(begin_block:end_block,begin_block:end_block) = cc{i}(begin_block:end_block,begin_block:end_block) + block;
          if precision
              cc{i} = inv(cc{i}) + eye(dim);
              %cc{i} = inv(cc{i});
          end
    end
end