function [X,V] = gestaltAncestralSample(ge,g,z,precision,varargin)
    parser = inputParser;
    addParameter(parser,'positive',false,@islogical);
    parse(parser,varargin{:});
    params = parser.Results;
    
    if ~precision
        Cv = componentSum(g,ge.cc);
    else
        Cv = inv(componentSum(g,ge.pc));
    end
    
    V = mvnrnd(zeros(ge.B,ge.Dv),Cv);
    if params.positive
        V = abs(V);        
    end

    means = reshape(V,ge.B,ge.Dv);
    means = z * means * ge.A';
    X = mvnrnd(means,ge.obsVar*eye(ge.Dx));
end