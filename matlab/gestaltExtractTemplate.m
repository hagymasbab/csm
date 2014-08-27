function template = gestaltExtractTemplate(ge,compNum)
    cc = ge.cc{compNum};
    thresh = (mean(cc(:)) + max(cc(:))) / 2 -1;
    template = zeros(ge.Dv,1);
    for i=1:ge.Dv
        for j=i+1:ge.Dv
            if cc(i,j) > thresh
                template(i,1) = 1;
                template(j,1) = 1;
            end
        end
    end
end