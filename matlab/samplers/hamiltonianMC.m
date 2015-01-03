function [samples,rr] = hamiltonianMC(init,logpdf,grad,nSamp,lfSteps,stepSize,varargin)
    % grad is the gradient of the logpdf, not of E, so it takes a negative
    % sign relative to Neal 2010.
    parser = inputParser;
    addParameter(parser,'verbose',0,@isnumeric);
    addParameter(parser,'plot',false,@islogical);
    addParameter(parser,'bounds',[]);
    parse(parser,varargin{:});
    v = parser.Results.verbose;
    pl = parser.Results.plot;
    bounds = parser.Results.bounds;
    colors = hsv(nSamp*10);
    
    E_grad = @(q) -grad(q);
    K = @(p) sum(p.^2);
    
    dim = size(init,1);
    samples = zeros(nSamp,dim);

    rr = 0;
    q_act = init;
    lp_act = logpdf(q_act);
    
    i = 1;
    while i<=nSamp
        % sample momentum
        p_act = randn(dim,1);
        
        % randomise stepsize to avoid (near) periodic trajectories
        act_step = (0.8 + rand()*0.4)*stepSize;
        % leapfrog integration of Hamiltonian dynamics
        [p,q,ptraj,traj] = leapfrog(p_act,q_act,E_grad,act_step,lfSteps,bounds);
        %q'
        %pause
        lp_next = logpdf(q);

        if pl
            subplot(2,1,1);
            plot(traj(1,:),traj(2,:),'color',colors(i,:))
            subplot(2,1,2);
            plot(sum(ptraj,1).^2);
            pause
        end
        
        % accept or reject
        m_act = sum(p_act.^2)/2;
        m_next = sum(p.^2)/2;
        H_act = -lp_act + m_act;
        H_next = -lp_next + m_next;
        a = rand(); 
        if v==2
            %fprintf('lp1 %.2e lp2 %.2e difflp %.2e m1 %.2e m2 %.2e diffm %.2e, diff %.2e ap. %.4f rand %.4f',lp_act,lp_next,lp_act-lp_next,m_act,m_next,m_act-m_next,lp_act+m_act-lp_next-m_next,exp(H_act-H_next),a);
            %fprintf('lp1 %.2f lp2 %.2f difflp %.2f m1 %.2f m2 %.2f diffm %.2f, diff %.2f ap. %.4f rand %.4f',lp_act,lp_next,lp_act-lp_next,m_act,m_next,m_act-m_next,lp_act+m_act-lp_next-m_next,exp(H_act-H_next),a);
            fprintf('lp1 %.2f lp2 %.2f ap. %.4f rand %.4f',lp_act,lp_next,exp(H_act-H_next),a);
        end
        if a < exp(H_act - H_next)
            samples(i,:) = q';
            q_act = q;
            lp_act = lp_next;
            
            if v==1
                printCounter(i,'stringVal','Sample','maxVal',nSamp,'newLine',true);
            elseif v==2
                fprintf(' accepted \n');
            end
            
            i = i+1;
        else
            if v==2
                fprintf(' rejected \n');
            end
            rr = rr+1;
        end
    end
end
        