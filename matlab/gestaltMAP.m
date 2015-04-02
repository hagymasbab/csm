function [V,G,Z] = gestaltMAP(ge,X)
    if ge.B ~= 1
        error('not implemented for B>1');
    end
    N = size(X,1);
    learning_rate_g = 0.0001;
    learning_rate_v = 0.01;
    learning_rate_z = 0.01;
    minval = 0.01;
    
    % initialise latent variables
    V = zeros(N,ge.Dv);
    V = (ge.A'*X')';
    G = 0.1*ones(N,ge.k);
    Z = ones(N,1);    
                   
    for i=1:N
        % gradient ascent
        convergence = false;
        while ~convergence
            %fprintf('.');
            % TODO reshape if needed
            grad = gestaltFullLogPosteriorGrad(ge,X(i,:),V(i,:),G(i,:)',Z(i,1),[]);
            %grad'
            if max(grad) > 1e-2
                grad_G = grad(1:ge.k,1);
                grad_V = grad(ge.k+1:end-1,1);
                grad_Z = grad(end,1);

                V(i,:) = V(i,:) + learning_rate_v * grad_V';
                G(i,:) = max(G(i,:) + learning_rate_g * grad_G',minval*ones(1,ge.k));
                Z(i,:) = max(Z(i,:) + learning_rate_z * grad_Z',minval);                
                %V(i,:)
                G(i,:)
                Z(i,:)
                pause
            else
                convergence = true;
            end
        end                
    end
end