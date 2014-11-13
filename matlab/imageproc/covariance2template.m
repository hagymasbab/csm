function template = covariance2template(cc)
    %cc = ge.cc{compNum};
    Dv = size(cc,1);
    cc = cc .* (ones(Dv)-eye(Dv));
    %thresh = mean(cc(:));
    thresh = (mean(abs(cc(:))) + max(abs(cc(:)))) / 2;
    %thresh = max(abs(cc(:))) / 2;
    %thresh = 0.01;
    template = zeros(Dv,1);
    for i=1:Dv
        for j=i+1:Dv
            if cc(i,j) > thresh
            %if abs(cc(i,j)) > thresh                
                template(i,1) = 1;
                template(j,1) = 1;
            end
        end
    end
end