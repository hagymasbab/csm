function [samples,rr] = gestaltSamplePosterior(ge,nSamp,varargin)
    parser = inputParser;
    addParamValue(parser,'verbose',0,@isnumeric);
    addParamValue(parser,'plot',0,@isnumeric);
    addParamValue(parser,'dnum',ge.N,@isnumeric);
    addParamValue(parser,'dshift',0,@isnumeric);
    addParamValue(parser,'precision',false,@islogical);
    addParamValue(parser,'savename','samples');
    parse(parser,varargin{:});
    v = parser.Results.verbose;
    pl = parser.Results.plot;
    precision = parser.Results.precision;   
    dnum = parser.Results.dnum;
    dshift = parser.Results.dshift + 1;
    % TODO test whether these values are valid
    fname = strcat(parser.Results.savename,'.mat');
    
    rr = 0;
    samples = zeros(dnum,nSamp,ge.k+(ge.B*ge.Dv));
    if v == 0
        fprintf('Observation %d/',dnum);
    end
    for n=dshift:dshift+dnum-1
        if v==0
            printCounter(n);
        end
        
        [s,rr_act] = gestaltGibbs(ge,n,nSamp,'slice',0.05,'verbose',v,'plot',pl,'burnin',0,'thin',1,'precision',precision);
        samples(n,:,:) = s;
        rr = rr + rr_act;
    end
    rr = rr/dnum;
    if v==0
        fprintf('\n');
    end
    % save samples
    save(fname,'samples');
end