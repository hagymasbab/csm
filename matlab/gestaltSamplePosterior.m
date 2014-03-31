function [samples,rr] = gestaltSamplePosterior(ge,nSamp,varargin)
    parser = inputParser;
    addParamValue(parser,'verbose',0,@isnumeric);
    addParamValue(parser,'plot',false,@islogical);
    addParamValue(parser,'dnum',ge.N,@isnumeric);
    addParamValue(parser,'dshift',0,@isnumeric);
    addParamValue(parser,'savename','samples');
    parse(parser,varargin{:});
    v = parser.Results.verbose;
    pl = parser.Results.plot;
    dnum = parser.Results.dnum;
    dshift = parser.Results.dshift + 1;
    % TODO test whether these values are valid
    fname = strcat(parser.Results.dnum,'.mat');

    init = [0.5*ones(ge.k,1); 0*ones(ge.Dv,1)];
    
    if pl
        clf;
        subplot(1,2,1);
        hold on;
        subplot(1,2,2);
        hold on;
    end
    
    rr = 0;
    samples = zeros(dnum,nSamp,ge.k+ge.Dv);
    if v == 0
        fprintf('Observation %d/',dnum);
    end
    for n=dshift:dshift+dnum-1
        if v==0
            printCounter(n);
        end
        % gestaltPlotPosterior(ge,ge.X(n,:)');
        %logpdf = @(gv) gestaltUnnormLogPost(gv(1:ge.k,:),gv(ge.k+1:size(gv),:),ge.X(n,:)',ge);
        %grad = @(gv) gestaltLogPostGradient(gv(1:ge.k,:),gv(ge.k+1:size(gv),:),ge.X(n,:)',ge,true);
        %gradObj = GestaltLogPostGrad(ge);
        %grad = @(gv) gradObj.statefulGrad(gv(1:ge.k,:),gv(ge.k+1:size(gv),:),ge.X(n,:)',ge);
        %coeff = symmetricDirichlet(ge.sparsity,ge.k,1)';
        %Cv = componentSum(coeff,ge.cc);
        %[s,rr_act] = combinedMC(init,1:ge.k,logpdf,grad,nSamp,Cv,0.05,'summedToOne',true,'plot',pl,'verbose',v);
        %[s,rr_act] = combinedMC(init,[],logpdf,grad,nSamp,Cv,0.05,'plot',pl,'verbose',v,'bounds',[1 0 1;2 0 1]);
        
        [s,rr_act] = gestaltGibbs(ge,n,nSamp,'slice',0.01,'verbose',v);
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