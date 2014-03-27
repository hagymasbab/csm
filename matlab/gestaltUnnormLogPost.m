function lp = gestaltUnnormLogPost(g,v,x,ge)
    if (sum(g<0) > 0) || (sum(g>1) > 0) 
        lp=-Inf;
%     elseif sum(g(1:ge.k-1,1)) == (1-g(ge.k,1)) 
%         lp=-Inf;
%         g
%         'eeeeeeeeeeeeeeeeeeee'
%         sum(g(1:ge.k-1,:))
%         1-g(ge.k,:)
%         %fprintf('proposal of g is negative or larger than 1');
    else
        %iCx = inv(ge.obsVar*eye(ge.Dx));
        Cx = ge.obsVar*eye(ge.Dx);
        Cv = componentSum(g,ge.cc);
        xAv = x-ge.A*v;
        prior = (ge.sparsity - 1) * sum(log(g));
        %lp = (-1/2) * (  xAv' * (Cx \  xAv) + logdet(Cv,'chol') + v'*(Cv \ v) ) + prior;
        lp = (-1/2) * (  xAv' * (Cx \  xAv) + log(abs(det(Cv))) + v'*(Cv \ v) ) + prior;
    end
end