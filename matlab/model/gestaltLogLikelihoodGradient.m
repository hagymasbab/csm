function grad = gestaltLogLikelihoodGradient(ge,L,data,varargin)
    parser = inputParser;   
    addParameter(parser,'verbose',0,@isnumeric);  
    addParameter(parser,'randseed','leave');      
    parse(parser,varargin{:});        
    params = parser.Results;  
    
    setrandseed(params.randseed);
    N = size(data,1);
    if ge.B ~= 1
        error('not implemented');
    end
    data = reshape(data,N,ge.Dx);
    
    V = gestaltSamplePriorV(ge,L,[]);
    Z = gamrnd(ge.z_shape,ge.z_scale,[L 1]);
    
    pA = pinv(ge.A);
    
end