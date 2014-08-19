function [samples,rr] = metropolisHastings(init,logpdf,propCov,nSamp,burnIn,thin,varargin)
    p = inputParser;
    addParamValue(p,'verbose',0,@isnumeric);
    addParamValue(p,'adapt',false,@islogical);
    addParamValue(p,'adaptStep',1.01,@isnumeric);
    addParamValue(p,'adaptSmooth',1,@isnumeric);
    p.KeepUnmatched = true;
    parse(p,varargin{:});
    v = p.Results.verbose;

    % the proposal density is a Gaussian
    N = nSamp * (thin + 1) + burnIn;
    d = size(init,1);
    allsamples = zeros(N,d); 
    
    act = init;
    lp_act = logpdf(act);
    if(v==2)
        fprintf('Initial log-prob: %e\n',lp_act);
    elseif v==1
        fprintf('Sample %d/', N);
    end
    reject = 0;
    propScale = 1;
    for i=1:N
        % draw a new sample from the proposal
        next = (mvnrnd(act,propScale*propCov))';
        lp_next = logpdf(next);
        a = rand();
        if(v==2)
            fprintf('Proposed log-prob: %e, acceptance prob %f, random num %f',lp_next, min(1,exp(lp_next) / exp(lp_act)),a);
        elseif v==1
            printCounter(i);
        end
        accept = lp_next > lp_act || log(a) < lp_next - lp_act;
        act_rej = 0;
        while ~accept
            if p.Results.adapt && act_rej > 5
                propScale = propScale / p.Results.adaptStep;
            end
            if v==2
                fprintf(' rejected\n');
            end
            act_rej = act_rej + 1;
            next = (mvnrnd(act,propCov))';
            lp_next = logpdf(next);
            a = rand();
            if(v==2)
                fprintf('Proposed log-prob: %e, acceptance prob %f, random num %f',lp_next, min(1,exp(lp_next) / exp(lp_act)),a);
            end
            accept = lp_next > lp_act || log(a) < lp_next - lp_act;
        end       
        if v==2
            fprintf(' accepted\n');
        end
        if p.Results.adapt
            % acceptance rate should be around 0.234
            if act_rej < 4
                % increase proposal variance
                propScale = propScale * p.Results.adaptStep;
            end
            if v==2
                fprintf('num of rejections: %d propsal scaled to %f\n',act_rej,propScale);
            end
        end
        reject = reject + act_rej;
        act = next;
        lp_act = lp_next;
        allsamples(i,:) = act';
    end
    rr = reject/N;
    if v==1
        fprintf(' rejections per acceptance on average %f\n',rr);
    end
    % burn-in and thin
    samples = allsamples;
    samples = samples(burnIn+1:size(samples,1),:);
    samples = samples(1:(thin+1):size(samples,1),:);

    