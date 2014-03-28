function [s,rr] = gestaltGibbs(ge,xind,nSamp,metsteps,metVar,varargin)
    parser = inputParser;
    addParamValue(parser,'verbose',0,@isnumeric);
    parse(parser,varargin{:});
    verb = parser.Results.verbose;    
    
    s = zeros(nSamp,ge.k + ge.Dv);
    rr = 0;
    g = 0.5 * ones(ge.k,1);
    v = zeros(ge.Dv,1);
    
    if verb==1
        fprintf('Sample %d/',nSamp);
    end
    for i=1:nSamp
        if verb==1
            printCounter(i);
        end
        
        % generate a direct sample from the conditional posterior over v
        v = gestaltPostVRnd(ge,xind,g);
        
        % Metropolis-Hastings scheme to sample from the conditional
        % posterior over g
        ms = 0;
        lp_act = gestaltLogPostG(g,v,ge);
        while ms < metsteps
            % propose from a unit Gaussian of dimension K-1
            g_part = mvnrnd(g(1:ge.k-1,1)',metVar*eye(ge.k-1))';
            % the last element is determined by the rest
            g_next = [g_part; 1-sum(g_part)];
            a = rand();
            lp_next = gestaltLogPostG(g_next,v,ge);
            if verb==2
                fprintf('%f %f %f %f ',lp_act,lp_next,exp(lp_next - lp_act),a);
                pause
            end
            
            if a < exp(lp_next - lp_act)
                % accept the sample
                g = g_next;
                lp_act = lp_next;
                ms = ms + 1;
                if verb==2
                    fprintf('accept %d\n',ms);
                end
            else
                rr = rr + 1;
                if verb==2
                    fprintf('reject\n');
                end
            end
        end
        
        %v = gestaltPostVRnd(ge,xind,g);
        
        % store the combined sample
        s(i,:) = [g' v'];
    end
    if verb==1
        fprintf('\n');
    end
    
    % calculate the rejection rate
    rr = rr / (rr + nSamp);
end