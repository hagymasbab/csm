function [samples,rr] = combinedMC(init,met_dims,logpdf,grad,nSamp,covariance,propVar,varargin)
    parser = inputParser;
    addParamValue(parser,'verbose',0,@isnumeric);
    addParamValue(parser,'plot',false,@islogical);
    addParamValue(parser,'summedToOne',false,@islogical);
    addParamValue(parser,'transformMomentum',false,@islogical);
    addParamValue(parser,'bounds',[]);
    parse(parser,varargin{:});
    v = parser.Results.verbose;
    pl = parser.Results.plot;
    bounds = parser.Results.bounds;     
    summed = parser.Results.summedToOne;     
    tm = parser.Results.transformMomentum;

    met = zeros(size(init));
    met(met_dims) = 1;
    dim = size(init,1);
    d_met = sum(met);
    d_ham = dim - d_met;
    
    % set step size and number of steps for leapfrog integration
    if tm
        % TODO
        mCov = covariance;
        iCov = inv(covariance);
    else
        mCov = eye(d_ham);
        iCov = eye(d_ham);
    end
    [~,S,~] = svd(covariance);
    leapsteps = ceil(max(diag(S)) / min(diag(S)));
    % TEST
    leapsteps = min(leapsteps,10000);
    stepsize = min(0.9*min(diag(S)),1);
    leapsteps = 100;
    stepsize = 0.02;
    if v>0
        fprintf('step size %f step num %d\n',stepsize,leapsteps);
    end

    propCov = propVar * eye(d_met-summed);
    % best to use a gradient function that does not do any computation for
    % the Met dimensions
    %E_grad = @(q) partialNegGrad(q,~met,grad);
    E_grad = @(q) -grad(q).*(~met);

    q_act = init;
    E_act = -logpdf(q_act);
    
    if v==1
        fprintf('Sample %d/',nSamp);
    end
    samples = zeros(nSamp,dim);
    rr = 0;
    color = 1;
    plothandle = 0;

    for i=1:nSamp
        met_act = q_act(met==1);
        
        % propose from a Gaussian for the Metropolis update until it could
        % be accepted without changing the other dimensions
        mtraj = [];
        if d_met>0 
            accepted = false;
            if v==2
                fprintf('Metropolis proposal ');
            end
            steps = 1;
            while ~accepted
                met_next = (mvnrnd(met_act(1:d_met-summed),propCov))';
                if summed
                    met_next = [met_next;1-sum(met_next)];
                end
                if v==2
                    printCounter(steps);
                end
                temporary_q = q_act;
                temporary_q(met==1) = met_next;
                met_next_E = -logpdf(temporary_q);
                b = rand();
                if b < exp(E_act - met_next_E)
                    accepted = true;
                    q_act = temporary_q;
                    E_act = met_next_E;
                    ind = size(mtraj)+1;
                    mtraj(:,ind:ind+1) = [met_act met_next];
                else
                    rr = rr + 1;
                    steps = steps + 1;
                end         
            end
            if v==2
                fprintf('\n');
            end
            if pl    
                subplot(1,2,1);
                if v==2
                    plot(mtraj(1,:),mtraj(2,:));
                end
                scatter(met_act(1,:),met_act(2,:));
            end
        end               
        
        % propose from the Hamiltonian dynamics for the rest
        ham_act = q_act(~met==1);
        if d_ham > 0
            accepted = false;
            while ~accepted
                %p_ham = mvnrnd(zeros(1,d_ham),eye(d_ham))';
                p_ham = mvnrnd(zeros(1,d_ham),iCov)';
                p_act = zeros(dim,1);
                p_act(~met==1) = p_ham;
                % randomise stepsize to avoid (near) periodic trajectories
                act_step = (0.8 + rand()*0.4)*stepsize;
                [p_next,q_next,~,qtraj] = leapfrog(p_act,q_act,E_grad,act_step,leapsteps,bounds);
                % assemble the full position vectors
                %q_next = q_act;
                %q_next(~met==1) = ham_next;

                % accept or reject according to the HMC rule using zero momenta
                % for the Metropolis dimensions
                E_next = -logpdf(q_next);
                %K_act = sum(p_act.^2)/2;
                %K_next = sum(p_next.^2)/2;
                K_act = (p_act(met==0)' * mCov * p_act(met==0)) / 2;
                K_next = (p_next(met==0)' * mCov * p_next(met==0)) / 2;
                a = rand();
                if v==2            
                    fprintf('E1 %.2f E2 %.2f diffE %.2f K1 %.2f K2 %.2f diffK %.2f, diff %.2f ap. %.4f rand %.4f',E_act,E_next,E_act-E_next,K_act,K_next,K_act-K_next,E_act+K_act-E_next-K_next,exp(E_act+K_act-E_next-K_next),a);
                end
                % plot
                if pl
                    subplot(1,2,2);                    
                    if v==2
                        plothandle = plot(qtraj(d_met+1,:),qtraj(d_met+2,:));
                        pause
                    end
                end
                if a < exp(E_act + K_act - E_next - K_next)
                    q_act = q_next;
                    E_act = -logpdf(q_act);
                    samples(i,:) = q_act';

                    if v==1
                        printCounter(i);
                    elseif v==2
                        fprintf(' accepted \n');
                    end
                    accepted = true;
                    if pl
                        scatter(q_act(d_met+1,:),q_act(d_met+2,:));
                    end
                else
                    rr = rr + 1;
                    if v==2
                        fprintf(' rejected \n');
                    end
                    
                end
                if pl && v==2
                    delete(plothandle);
                end
            end % while ~accepted
        end % if d_ham > 0            
    end % for i=1:nSamp
    
    if v==1
        fprintf('\n');
    end
    
    rr = rr/(rr+nSamp);
    
    function dq = partialNegGrad(q,dims,grad)
        state = zeros(size(dims));
        state(dims==1) = q;
        dq = -grad(state);
        dq = dq(dims==1);
    end
end
