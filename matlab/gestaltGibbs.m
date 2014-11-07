function [s,rr,zsamp] = gestaltGibbs(ge,xind,nSamp,varargin)
    parser = inputParser;
    addParamValue(parser,'verbose',0,@isnumeric);
    addParamValue(parser,'stepsize',0.05,@isnumeric);
    addParamValue(parser,'burnin',0,@isnumeric);
    addParamValue(parser,'thin',1,@isnumeric);
    addParamValue(parser,'plotG',false,@islogical);
    addParamValue(parser,'plotZ',false,@islogical);
    addParamValue(parser,'sampleRetry',10,@isnumeric);
    addParamValue(parser,'precision',false,@islogical);
    addParamValue(parser,'contrast',false,@islogical);
    addParamValue(parser,'initG',[]);
    addParamValue(parser,'priorG','gamma');
    addParamValue(parser,'gSampler','gibbs-slice');
    addParamValue(parser,'repeatCycle',1,@isnumeric);
    parse(parser,varargin{:});
    params = parser.Results;

    N = nSamp*params.thin + params.burnin;        
    
    s = zeros(N,ge.k + ge.B*ge.Dv);
    zsamp = zeros(N,1);
    rr = 0;

    if isempty(params.initG)
        % sample g from the prior distribution
        try
            g = gestaltSamplePriorG(ge,params.priorG,'sampleRetry',params.sampleRetry);                   
        catch err
            rethrow(err);
        end      
    else
        g = params.initG;
    end

    V = zeros(ge.B,ge.Dv); % unused if we sample the conditional over v first
    z = 1; % remains unchanged if we do not use a contrast variable
    % TODO we might sample z from its prior
        
    for i=1:N
        if params.verbose==1
            printCounter(i,'stringVal','Sample','maxVal',N,'newLine',false);
        end
        
        % generate a direct sample from the conditional posterior over v        
        V = gestaltPostVRnd(ge,xind,g,z,params.precision);
        
%         if params.plot > 0
%             clf;
%             gestaltPlotCondPostG(ge,V,params.priorG,0.1);
%             hold on;
%             plot(ge.G(xind,1),0,'go');
%             pause
%         end
        
        % slice sampling for g
        logpdf = @(g) gestaltLogPostG(g,V,ge,params.priorG,params.precision); 
        
        valid = false;
        tries = 0;
        while ~valid
            if tries > params.sampleRetry
                if params.verbose == 1
                    delPrint(i);
                end
                throw(MException('Gestalt:Gibbs:TooManyTries','Number of tries to sample a valid g vector from the conditional posterior exceeded %d at sampling step %d',params.sampleRetry,i));
            end
            tries = tries + 1;
            
            try
                if strcmp(params.gSampler,'slice')
                    if strcmp(params.priorG,'dirichlet')
                        [g_part,rr_act] = sliceSample(g(1:ge.k-1,1),logpdf,params.stepsize,'plot',params.plot>1,'limits',[0,1]);
                        g_temp = [g_part; 1-sum(g_part)];

                    else
                        [g_temp,rr_act] = sliceSample(g,logpdf,params.stepsize,'plot',params.plot>1);
                    end
                elseif strcmp(params.gSampler,'mh')
                        [g_temp,rr_act] = metropolisHastings(g,logpdf,0.001*eye(ge.k),1,0,0,'verbose',0);
                        g_temp = g_temp';
                elseif strcmp(params.gSampler,'gibbs-slice')
                    if strcmp(params.priorG,'dirichlet')
                        throw(MException('Gestalt:Gibbs:InvalidParameterCombination','Gibbs-slice sampling is only implemented for independent (e.g. Gamma) priors.'));
                    else
                        rr_act = 0;
                        g_temp = g;
                        for cycle = 1:params.repeatCycle
                            for j = 1:ge.k
                                if params.plotG
                                    clf;
                                    gestaltPlotConditional(g,j,V,ge,params.priorG,0.1);
                                    hold on;
                                    scatter(ge.G(xind,j),0,140,'go','LineWidth',3);
                                    pause
                                end
                                condlogpdf = @(gi) gestaltLogCondPostG(gi,g_temp,j,V,ge,params.priorG,params.precision); 
                                [g_temp(j,1),rr_part] = sliceSample(g_temp(j,1),condlogpdf,params.stepsize,'plot',params.plotG);

                                rr_act = rr_act + rr_part;
                            end
                        end
                    end
                end
            catch err
                if strcmp(err.identifier,'Gestalt:SliceSample:TooManyTries')
                    fprintf('too many retries in slice sampling at datum %d sample %d try %d\n',xind,i,tries);
                    continue
                else
                    error(err.message);
                end
            end
                               
            valid = gestaltCheckG(g_temp,ge,params.precision);            
            if valid
                g = g_temp;
            else
                %fprintf('invalid');
            end
        end
        rr = rr + rr_act;
        
        % slice sampling for z
        if params.contrast            
            zlogpdf = @(z) gestaltLogPostZ(z,xind,V,ge); 
            valid = false;
            tries = 0;
            while ~valid
                if tries > params.sampleRetry    
                    throw(MException('Gestalt:Gibbs:TooManyTries','Number of tries to sample a valid z from the conditional posterior exceeded %d at sampling step %d',params.sampleRetry,i));
                end
                
                valid = true;                
                if params.plotZ
                    clf;                    
                    gestaltPlotZCond(ge,xind,V);
                    title(sprintf('N=%d L=%d try %d',xind,i,tries+1));
                    hold on;
                    %scatter(ge.Z(xind,1),0,140,'go','LineWidth',3);
                    pause
                end
                try
                    [z,rr_act] = sliceSample(z,zlogpdf,params.stepsize,'plot',params.plotZ,'limits',[0,Inf]);
                catch err
                    if strcmp(err.identifier,'Gestalt:SliceSample:TooManyTries')
                        valid = false;
                    else
                        throw(err);
                    end
                end
                tries = tries + 1;                       
            end
            rr = rr + rr_act;
        end
        
        % store the combined sample
        vlong = reshape(V,1,ge.B*ge.Dv);
        s(i,:) = [g' vlong];
        zsamp(i,1) = z;
    end
    
    % calculate the rejection rate
    rr = rr / (rr + N);
    
    % discard burn-in stage
    if params.burnin > 0
        s = s(params.burnin+1:N,:);
    end
    if params.thin > 1
        indices = 1:params.thin:nSamp*params.thin;
        s = s(indices,:);
    end
end