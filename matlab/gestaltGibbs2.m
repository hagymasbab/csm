function [vsamp,gsamp,zsamp,rr] = gestaltGibbs2(ge,xind,nSamp,varargin)    
    if ge.B > 1
        exit('Gibbs2 is not implemented for B>1');
    end

    parser = inputParser;
    addParameter(parser,'burnin',0,@isnumeric);
    addParameter(parser,'thin',1,@isnumeric);
    addParameter(parser,'stepsize',0.5,@isnumeric);
    addParameter(parser,'verbose',0,@isnumeric);
    parse(parser,varargin{:});
    params = parser.Results;

    N = nSamp*params.thin + params.burnin;   
    
    vsamp = zeros(nSamp,ge.Dv);
    gsamp = zeros(nSamp,ge.k);
    zsamp = zeros(nSamp,1);
    
    % initialisation
    
    z = 1;
    g = gestaltSamplePriorG(ge,1)';    
    Cv = componentSum(g,ge.cc);
    
    rr = 0;
    
    ge.sparseComponents = length(find(componentSum(1,ge.cc))) < ge.Dv^2 / 2;       
    
    for i=1:N
        if params.verbose==1
            printCounter(i,'stringVal','Sample','maxVal',N);
        end
                
        [~,e] = chol(Cv);
        if e ~= 0
            Cv = nearestSPD(Cv);
        end
        
        % sample V
        V = gestaltPostVRnd(ge,xind,Cv,z,false);
        
        % sample Z
        zlogpdf = @(z) gestaltLogPostZ(z,xind,V,ge);
        [z,rr_act] = sliceSample(z,zlogpdf,params.stepsize,'limits',[0,Inf]);        
        
        % sample G
        g_temp = g;
        for j = 1:ge.k
            prev_g = g_temp;
            condlogpdf = @(gi) gestaltLogCondPostG(gi,g_temp,j,V,ge,ge.prior,false,Cv,ge.cc); 
            [g_temp(j,1),rr_part] = sliceSample(g_temp(j,1),condlogpdf,params.stepsize,'limits',[0,Inf]);

            rr_act = rr_act + rr_part;
            actdiff = prev_g(j,1) - g_temp(j,1);
            Cv = Cv - actdiff * ge.cc{j};
        end
        g = g_temp;
        
        rr = rr + rr_act;
        
        if i > params.burnin && rem((i-params.burnin),params.thin) == 0
            vsamp(i,:) = V;
            gsamp(i,:) = g';
            zsamp(i,1) = z;
        end
    end
    
    rr = rr / (rr + N);
    
end
            