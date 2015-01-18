function printProgress(n,name)
    if nargin ==2
        fprintf('\n');
        for i = 1:n
            if length(name)>=i && i <= n-1
                fprintf(name(i));
            else
                fprintf('.');
            end
        end
        fprintf(']\n\n');
    else
        fprintf('\b|\n');
    end
end