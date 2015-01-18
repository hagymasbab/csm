function printProgress(n,name)
    fprintf('\n');
    for i = 1:n-1
        if length(name)>=i && i <= n-2
            fprintf(name(i));
        else
            fprintf('.');
        end
    end
    fprintf(']\n\n');
end