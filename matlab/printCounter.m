function printCounter(i)
    if i>1
          for j=0:log10(i-1);
              fprintf('\b'); % delete previous counter display
          end
    end
    fprintf('%d', i);
end