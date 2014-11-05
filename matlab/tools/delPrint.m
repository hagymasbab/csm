function delPrint(num)
    for b=1:9+2*(floor(log10(num))+1)
        fprintf('\b');
    end
end