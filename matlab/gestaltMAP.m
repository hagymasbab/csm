function [V,G,Z] = gestaltMAP(ge,X)
    close all;
%     if ge.B ~= 1
%         error('not implemented for B>1');
%     end
    N = size(X,1);
    learning_rate_g = 0.0001;
    learning_rate_v = 0.01;
    learning_rate_z = 0.001;
    minval = 0.01;
    
    % initialise latent variables
    %V = zeros(N,ge.Dv);
    %V = (ge.A'*X')';
    V = ge.V;
    %G = (1/ge.k)*ones(N,ge.k);
    G = rand(N,ge.k)
    %G = repmat([0.01 1],N,1);
    %Z = 1 * ones(N,1);                           
    Z = ge.Z;
    
    for i=1:N
        printCounter(i,'maxVal',N,'stringVal','Observation');        
        if ge.B==1
            act_x = X(i,:);
            act_v = V(i,:);
            norm_v = reshape(ge.V(i,1,:),1,ge.Dv);
            norm_v = norm_v / norm(norm_v);
        else
            act_x = X(i,:,:);
            act_v = V(i,:,:);
        end
        
        act_g = G(i,:)';
        act_z = Z(i,1);
        %lp = gestaltFullLogPosterior(ge,act_x,act_v,act_g,act_z,[]);
        % gradient ascent
        convergence = false;
        j = 1;
        while ~convergence
            %printCounter(j);
            j = j+1;
            %fprintf('.');
            % TODO reshape if needed
            if ge.B==1
                grad = gestaltFullLogPosteriorGrad(ge,act_x,act_v,act_g,act_z,[]);
                grad_G = grad(1:ge.k,1)
                grad_V = grad(ge.k+1:end-1,1);
                grad_Z = grad(end,1);
                
                pref_shift = 0.1;
                if max(abs(grad_G)) * learning_rate_g > pref_shift                    
                    learning_rate_g = pref_shift / max(abs(grad_G))
                    learning_rate_z = pref_shift / max(abs(grad_Z));
                    learning_rate_v = pref_shift / max(abs(grad_V));
                end
            end
            
%             gx = 0.01:0.02:1.1;
%             pg = zeros(length(gx));
%             for k=1:length(gx);printCounter(k,'maxVal',length(gx),'StringVal','e');for j=1:length(gx);pg(k,j) = gestaltFullLogPosterior(ge,reshape(act_x,ge.B,ge.Dx),reshape(act_v,ge.B,ge.Dv),[gx(k);gx(j)],act_z,[]);end;end;
%             pg (pg < max(pg(:))*0.99) = max(pg(:))*0.99;
%             cutoff = 1;
%             imagesc(gx(cutoff:end),gx(cutoff:end),pg(cutoff:end,cutoff:end))
            
            %max(grad)
            %min(grad)
            %grad'
            if max(abs(grad)) > 1e-1

%                 zx = 0.01:0.02:7;
%                 pz = zeros(length(zx),1);
%                 %for j=1:length(zx);pz(j) = gestaltFullLogPosterior(ge,act_x,ge.V(1,:,:),act_g,zx(j),[]);end;
%                 for j=1:length(zx);pz(j) = gestaltFullLogPosterior(ge,act_x,act_v,act_g,zx(j),[]);end;
%                 plot(zx,pz)                                
                
%                 gx = 0.01:0.1:2;
%                 gg1 = zeros(length(gx));
%                 gg2 = zeros(length(gx));
%                 for k=1:length(gx)
%                     printCounter(k,'maxVal',length(gx),'StringVal','e')
%                     for j=1:length(gx);
%                         gradient = gestaltFullLogPosteriorGrad(ge,act_x,act_v,[gx(k);gx(j)],act_z,[]);
%                         gg1(k,j) = gradient(1);
%                         gg2(k,j) = gradient(2);
%                     end
%                 end
%                 cutoff = 3;
%                 figure
%                 imagesc(gx(cutoff:end),gx(cutoff:end),gg1(cutoff:end,cutoff:end))
%                 title('1')
%                 figure
%                 imagesc(gx(cutoff:end),gx(cutoff:end),gg2(cutoff:end,cutoff:end))
%                 title('2')
                
                max_shift_g = 0.1;
                act_v = act_v + learning_rate_v * grad_V';
                act_g = max(act_g + min(learning_rate_g * grad_G,max_shift_g),minval*ones(ge.k,1))
                act_z = max(act_z + learning_rate_z * grad_Z,minval);
                lp = gestaltFullLogPosterior(ge,act_x,act_v,act_g,act_z,[]);                
                
%                 actnorm = reshape(act_v,1,ge.Dv);
%                 actnorm = actnorm/norm(actnorm);
%                 actnorm * norm_v';
%                 %viewImage(grad_V)
%                 fprintf('g %.2f v %.2f z %.2f lp %.2f\n',max(grad_G),max(grad_V),max(grad_Z),lp);
                pause
            else
                convergence = true;
            end
        end
        
        V(i,:) = act_v;
        G(i,:) = act_g';
        Z(i,1) = act_z;
    end
end