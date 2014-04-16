function [s,rr] = gestaltGibbs(ge,xind,nSamp,stepsize,varargin)
    parser = inputParser;
    addParamValue(parser,'verbose',0,@isnumeric);
    addParamValue(parser,'burnin',0,@isnumeric);
    addParamValue(parser,'thin',1,@isnumeric);
    addParamValue(parser,'plot',0,@isnumeric);
    addParamValue(parser,'precision',false,@islogical);
    parse(parser,varargin{:});
    verb = parser.Results.verbose;    
    pl = parser.Results.plot;    
    precision = parser.Results.precision;
    burn = parser.Results.burnin;
    thin = parser.Results.thin;
    N = nSamp*thin + burn;
    
    s = zeros(N,ge.k + ge.B*ge.Dv);
    rr = 0;
    %g = 0.5 * ones(ge.k,1);
    valid = false;
    while ~valid
        %fprintf('e');
        g = symmetricDirichlet(0.2,ge.k,1)';
        valid = checkG(g,ge,precision);
    end
    V = zeros(ge.B,ge.Dv); % unused if we sample the conditional over v first
        
    if verb==1
        fprintf('Sample %d/',N);
    end
    for i=1:N
        if verb==1
            printCounter(i);
        end
        
        % generate a direct sample from the conditional posterior over v        
        V = gestaltPostVRnd(ge,xind,g,precision);
        
        if pl > 0
            clf;
            gestaltPlotCondPostG(ge,V);
            hold on;
            plot(ge.G(xind,1),0,'ro');
            pause
        end
        
        % slice sampling for g
        if ~precision
            logpdf = @(g) gestaltLogPostG(g,V,ge); 
        else
            logpdf = @(g) gestaltLogPostGPrec(g,V,ge); 
        end
        valid = false;
        while ~valid
            [g_part,rr_act] = sliceSample(g(1:ge.k-1,1),logpdf,stepsize,'plot',pl>1);
            g = [g_part; 1-sum(g_part)];
            valid = checkG(g,ge,precision);
        end
        rr = rr + rr_act;

        % uncomment this and comment out the similar line in the beginning if
        % you want to reverse the order of sampling from the conditionals
        % if ~precision
        %     V = gestaltPostVRnd(ge,xind,g);
        % else
        %     V = gestaltPostVRndPrec(ge,xind,g);
        % end
        
        % store the combined sample
        vlong = reshape(V,1,ge.B*ge.Dv);
        s(i,:) = [g' vlong];
    end
%     if verb==1
%         fprintf('\n');
%     end
    
    % calculate the rejection rate
    rr = rr / (rr + N);
    
    % discard burn-in stage
    if burn > 0
        s = s(burn+1:N,:);
    end
    if thin > 1
        indices = 1:thin:nSamp*thin;
        s = s(indices,:);
    end
end

function good = checkG(g,ge,precision)
    if ~precision
        CvP = componentSum(g,ge.cc);
    else
        CvP = componentSum(g,ge.pc);
    end
    if rcond(CvP) < 1e-15
        good = false;
        return;
    end
    if ~precision
        postP = (1/ge.obsVar) * ge.AA + inv(CvP);                
    else
        postP = (1/ge.obsVar) * ge.AA + CvP;                
    end
    if rcond(postP) < 1e-15
        good = false;
        return;
    end
    postC = inv(postP);
    [~,err] = cholcov(postC);
    if det(CvP) > 0 && det(postC) > 0 && err == 0                
        good = true;                
    else
        good = false;
    end
end