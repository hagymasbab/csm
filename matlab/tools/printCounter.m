function printCounter(i,varargin)
    parser = inputParser;    
    addParamValue(parser,'maxVal',0,@isnumeric);
    addParamValue(parser,'stringVal','');
    addParamValue(parser,'newLine',true,@islogical);
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