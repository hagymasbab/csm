function cc = randomCovariances(k,dim)
    blocksize = floor(dim/k);
    block = 0.5 * ones(blocksize);
    for i=1:k
          base = rand(dim,dim);
          cc{i} = (1/dim)*(base'*base);
          cc{i} = cc{i} + eye(dim);
          begin_block = (i-1)*blocksize + 1;
          end_block = i*blocksize;
          cc{i}(begin_block:end_block,begin_block:end_block) = cc{i}(begin_block:end_block,begin_block:end_block) + block;
    end
end