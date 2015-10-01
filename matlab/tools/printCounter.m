function printCounter(i,varargin)
% printCounter - displays a simple counter within a for loop
% usage: for i=1:imax; printCounter(i,'maxVal',imax,'stringVal','Iteration','newLine',true); end
%
% Arguments:
%   i (integer) - the loop variable
%   Optional:
%     'maxVal' (integer) - see the coutner in a form maxVal/i, e.g. 100/23
%     'stringVal' (string) - display a string in front of the counter, e.g. Iter 100/23
%     'newLine' (logical) - whether to put the cursor in a new line at the end of the loop

    parser = inputParser;    
    addParameter(parser,'maxVal',0,@isnumeric);
    addParameter(parser,'stringVal','');
    addParameter(parser,'newLine',true,@islogical);
    parse(parser,varargin{:});
    params = parser.Results;
    
    if ~isempty(params.stringVal) && params.maxVal > 0 && i == 1
        fprintf('%s %d/',params.stringVal,params.maxVal);
    end
    
    if i>1
          for j=0:log10(i-1);
              fprintf('\b'); % delete previous counter display
          end
    end
    fprintf('%d', i);
    
    if i == params.maxVal
        if params.newLine
            fprintf('\n');
        else
            fprintf(' ');
        end        
    end
end