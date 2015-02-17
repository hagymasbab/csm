function [X,V] = gestaltAncestralSample(ge,g,z,varargin)
    parser = inputParser;
    addParameter(parser,'positive',false,@islogical);
    addParameter(parser,'precision',false,@islogical);
    addParameter(parser,'N',1,@isnumeric);
    parse(parser,varargin{:});
    params = parser.Results;
    
    if ~params.precision
        Cv = componentSum(g,ge.cc);
    else
        Cv = inv(componentSum(g,ge.pc));
    end
    
    nsamp = ge.B;
    if ge.B == 1
        nsamp = params.N;
    end
        
    
    V = mvnrnd(zeros(nsamp,ge.Dv),Cv);
    if params.positive
        V = abs(V);        
    end

    means = reshape(V,nsamp,ge.Dv);
    means = z * means * ge.A';
    X = mvnrnd(means,ge.obsVar*eye(ge.Dx));
end