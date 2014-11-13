function template = gestaltExtractTemplate(ge,compNum)
    cc = ge.cc{compNum};
    cc = cc .* (ones(ge.Dv)-eye(ge.Dv));
    %thresh = mean(cc(:));
    %thresh = (mean(cc(:)) + max(abs(cc(:)))) / 2;
    thresh = max(abs(cc(:))) / -0.1;
    %thresh = 0.01;
    template = zeros(ge.Dv,1);
    for i=1:ge.Dv
        for j=i+1:ge.Dv
            if abs(cc(i,j)) > thresh
                template(i,1) = 1;
                template(j,1) = 1;
            end
        end
    end
end