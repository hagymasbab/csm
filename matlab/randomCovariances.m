function cc = randomCovariances(k,dim)
    df = 1023;
    for i=1:k
          base = rand(dim,dim);
          cc{i} = (1/dim)*(base'*base);
          cc{i} = cc{i} + eye(dim);
    end
end