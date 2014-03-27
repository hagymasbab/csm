function cc = randomCovariances(k,dim)
    df = 1023;
    for i=1:k
          base = rand(dim,dim);
          cc{i} = (1/dim)*(base'*base);
%         fprintf('Component %d ',i);
%         invertible = false;
%         rejected = 0;
%         while ~invertible
%              printCounter(rejected);
%             cc{i} =  wishrnd(base,df)/df;
%             if rcond(cc{i})>1e-16
%                 invertible = true;
%             else
%                 rejected = rejected + 1;
%             end
%         end
%        fprintf('\n');
    end
end