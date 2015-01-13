function [vsamp,gsamp,zsamp,rr] = gestaltGibbs(ge,xind,nSamp,varargin)
    parser = inputParser;
    addParameter(parser,'verbose',0,@isnumeric);
    addParameter(parser,'stepsize',0.1,@isnumeric);
    addParameter(parser,'burnin',0,@isnumeric);
    addParameter(parser,'thin',1,@isnumeric);
    addParameter(parser,'plotG',false,@islogical);
    addParameter(parser,'plotZ',false,@islogical);
    addParameter(parser,'sampleRetry',10,@isnumeric);
    addParameter(parser,'precision',false,@islogical);
    addParameter(parser,'contrast',true,@islogical);
    addParameter(parser,'initG',[]);
    addParameter(parser,'gSampler','gibbs-slice');
    addParameter(parser,'zSampler','slice');
    addParameter(parser,'vSampler','direct');
    addParameter(parser,'repeatCycle',1,@isnumeric);
    addParameter(parser,'prestimSamples',0,@isnumeric);
    addParameter(parser,'poststimSamples',0,@isnumeric);
    addParameter(parser,'initZ',1,@isnumeric);
    parse(parser,varargin{:});
    params = parser.Results;

    N = nSamp*params.thin + params.burnin;        
    
    % TODO do something when there is burn-in and prestim at the same time
    
    s = zeros(N,ge.k + ge.B*ge.Dv);
    zsamp = zeros(N,1);
    rr = 0;

    if isempty(params.initG)
        % sample g from the prior distribution
        try
            g = gestaltSamplePriorG(ge,ge.prior,'sampleRetry',params.sampleRetry);                   
        catch err
            rethrow(err);
        end      
    else
        g = params.initG;
    end
    
    
    V = zeros(ge.B,ge.Dv); % unused if we sample the conditional over v first
    z = params.initZ; % remains unchanged if we do not use a contrast variable
    % TODO we might sample z from its prior
    
    nullStimulus = ge.obsVar * randn([1,ge.B,ge.Dx]);
    storedStimulus = ge.X(xind,:,:);
    ge.X(xind,:,:) = nullStimulus;
    
    switched = 0;
    for i=1:N
        if params.verbose==1
            printCounter(i,'stringVal','Sample','maxVal',N,'newLine',false);
        end
        
        if i > params.prestimSamples && switched == 0
            ge.X(xind,:,:) = storedStimulus;
            switched = 1;
        elseif i > N - params.poststimSamples && switched == 1
            ge.X(xind,:,:) = nullStimulus;
            switched = 2;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
        % SAMPLE V       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % generate a direct sample from the conditional posterior over v  
        if strcmp(params.vSampler,'direct')
            V = gestaltPostVRnd(ge,xind,g,z,params.precision);
        elseif strcmp(params.vSampler,'test')
            V = ge.V(xind,:,:);
        elseif strcmp(params.vSampler,'mh')
            if ge.B ~= 1
                throw(MException('Gestalt:Gibbs:NotImplemented','Gibbs-MH sampling for V is not implemented for B > 0'));
            end
            %act_v = reshape(V,1,ge.Dv);
            sAA = ((z*z)/ge.obsVar) * ge.AA;
            Cv = componentSum(g,ge.cc);
            cov = inv(sAA + inv(Cv));  
            ATx = ge.A' * reshape(ge.X(xind,1,:),ge.Dv,1);
            m = ((z/ge.obsVar) * cov * ATx)';
            vlogpdf = @(v) log( mvnpdf(v',m,cov) );
            propcov = 0.005 * cov;
            %propcov = 0.00005 * eye(ge.Dv);
            [v_temp,rr_act] = metropolisHastings(V',vlogpdf,propcov,1,0,0,'verbose',0);
            V = reshape(v_temp,1,ge.Dv);
            rr = rr + rr_act;
        end
        
%         if params.plot > 0
%             clf;
%             gestaltPlotCondPostG(ge,V,ge.prior,0.1);
%             hold on;
%             plot(ge.G(xind,1),0,'go');
%             pause
%         end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
        % SAMPLE G       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % slice sampling for g
        logpdf = @(g) gestaltLogPostG(g,V,ge,ge.prior,params.precision); 
        
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
                    if strcmp(ge.prior,'dirichlet')
                        [g_part,rr_act] = sliceSample(g(1:ge.k-1,1),logpdf,params.stepsize,'plot',params.plotG,'limits',[0,1]);
                        g_temp = [g_part; 1-sum(g_part)];

                    else
                        [g_temp,rr_act] = sliceSample(g,logpdf,params.stepsize,'plot',params.plotG,'limits',[0,Inf]);
                    end
                elseif strcmp(params.gSampler,'mh')
                        [g_temp,rr_act] = metropolisHastings(g,logpdf,0.001*eye(ge.k),1,0,0,'verbose',0);
                        g_temp = g_temp';
                elseif strcmp(params.gSampler,'gibbs-slice')
                    if strcmp(ge.prior,'dirichlet')
                        throw(MException('Gestalt:Gibbs:InvalidParameterCombination','Gibbs-slice sampling is only implemented for independent (e.g. Gamma) priors.'));
                    else
                        rr_act = 0;
                        g_temp = g;
                        for cycle = 1:params.repeatCycle
                            cv = componentSum(g,ge.cc);                                                        
                            for j = 1:ge.k
                                if params.plotG
                                    clf;
                                    gestaltPlotConditional(g,j,V,ge,ge.prior,0.1);
                                    hold on;
                                    scatter(ge.G(xind,j),0,140,'go','LineWidth',3);
                                    pause
                                end                                
                                
                                prev_g = g_temp;
                                condlogpdf = @(gi) gestaltLogCondPostG(gi,g_temp,j,V,ge,ge.prior,params.precision,cv,ge.cc); 
                                [g_temp(j,1),rr_part] = sliceSample(g_temp(j,1),condlogpdf,params.stepsize,'plot',params.plotG,'limits',[0,Inf]);

                                rr_act = rr_act + rr_part;
                                actdiff = prev_g(j,1) - g_temp(j,1);
                                cv = cv - actdiff * ge.cc{j};
                            end
                        end
                    end
                elseif strcmp(params.gSampler,'test')
                    g_temp = ge.G(xind,:)';
                    rr_act = 0;
                end
            catch err
                if strcmp(err.identifier,'Gestalt:SliceSample:TooManyTries')
                    fprintf('too many retries in sampling of G at datum %d sample %d try %d\n',xind,i,tries);
                    continue
                else
                    %err.stack
                    error(err.getReport());
                    
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
        % SAMPLE Z       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
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
                    if strcmp(params.zSampler,'slice')
                        [z,rr_act] = sliceSample(z,zlogpdf,params.stepsize,'plot',params.plotZ,'limits',[0,Inf]);
                    elseif strcmp(params.zSampler,'test')
                        z = ge.Z(xind,1);
                        rr_act = 0;
                    end
                catch err
                    if strcmp(err.identifier,'Gestalt:SliceSample:TooManyTries')
                        fprintf('too many retries in sampling of Z at datum %d sample %d try %d\n',xind,i,tries);
                        continue
                    else
                        error(err.message);
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
    
    % TODO eliminate this, here because of learning
    if nargout > 2
        gsamp = s(:,1:ge.k);
        vsamp = reshape(s(:,ge.k+1:end),[nSamp ge.B ge.Dv]);    
    else
        vsamp = s;
        gsamp = rr;
    end
end