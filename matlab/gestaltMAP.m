function [V,G,Z,loglike] = gestaltMAP(ge,fix_v,fix_z,fix_g,v_init,drawGradients,saveCourses)
%function [V,G,Z,delta,gcourse,zcourse,loglike,initV] = gestaltMAP(ge,fix_v,fix_z,v_init,drawGradients)
    close all;
     if ge.B ~= 1
         error('not implemented for B>1');
     end
    X = reshape(ge.X,ge.N,ge.Dx);
    N = size(X,1);
    
    if fix_v
        learning_rate_v = 0;
    else
        learning_rate_v = 0.01;
    end
    
    if strcmp(v_init,'data')
        V = (ge.A'*X')';
    elseif strcmp(v_init,'true')
        V = ge.V;
    elseif strcmp(v_init,'rand')
        V = randn(N,ge.Dv);
    else
        V = v_init;
    end
    
    initV = V;
    
    if fix_z
        learning_rate_z = 0;
        Z = ge.Z;
    else
        learning_rate_z = 0.001;
        %Z = 1 * ones(N,1);
        Z = gamrnd(2*ones(N,1),2);
    end
    
    if fix_g
        learning_rate_g = 0;
        G = ge.G;
    else
        learning_rate_g = 0.001;
        G = rand(N,ge.k);
    end
    minval = 0.01;
    
    if saveCourses
        vcourse = {};
        gcourse = {};
        zcourse = {};               
    end
    loglike = {};
    
    for i=1:N
        
        act_lr_g = learning_rate_g;
        act_lr_v = learning_rate_v;
        act_lr_z = learning_rate_z;
        delta_threshold = 1e-6;
        
        act_x = X(i,:);
        act_v = V(i,:);
%         norm_v = reshape(ge.V(i,1,:),1,ge.Dv);
%         norm_v = norm_v / norm(norm_v);
        
        act_g = G(i,:)';
        act_z = Z(i,1);

        % gradient ascent
        convergence = false;
        counter = 1;
        
        if saveCourses
            actvcourse = [];
            actgcourse = [];
            actzcourse = [];
            actloglike = [];
        end
        while ~convergence       
            printCounter(counter)
            counter = counter+1; 

            grad = gestaltFullLogPosteriorGrad(ge,act_x,act_v,act_g,act_z,[]);
            grad_G = grad(1:ge.k,1);
            grad_V = grad(ge.k+1:end-1,1);
            grad_Z = grad(end,1);
%                 
%                 pref_shift = 0.1;
%                 if max(abs(grad_G)) * learning_rate_g > pref_shift                    
%                     learning_rate_g = pref_shift / max(abs(grad_G))
%                     learning_rate_z = pref_shift / max(abs(grad_Z))
%                     learning_rate_v = pref_shift / max(abs(grad_V))
%                 end
                
            if rem(counter-2,drawGradients) == 0
                imgmax = 2;
                %imgstep = 0.04;                    
                imgres = 30;
                probcutoff = 1;
                cutoff = 1;

                %gx = 0.01:imgstep:imgmax;
                gx = linspace(minval,imgmax,imgres);
                if max(act_g) > imgmax
                    gx = gx + max(act_g) - imgmax/2;
                end
                gg1 = zeros(length(gx));
                gg2 = zeros(length(gx));
                for k=1:length(gx)
                    printCounter(k,'maxVal',length(gx),'StringVal','e')
                    for j=1:length(gx);
                        gradient = gestaltFullLogPosteriorGrad(ge,act_x,act_v,[gx(k);gx(j)],act_z,[]);
                        gg1(k,j) = gradient(1);
                        gg2(k,j) = gradient(2);
                    end
                end           

                if ~exist('fhg1') 
                    fhg1 = figure;
                else 
                    figure(fhg1);
                end

                subplot(1,4,1);
                pg = zeros(length(gx));
                for k=1:length(gx)
                    printCounter(k,'maxVal',length(gx),'StringVal','e');
                    for j=1:length(gx);
                        pg(k,j) = gestaltFullLogPosterior(ge,reshape(act_x,ge.B,ge.Dx),reshape(act_v,ge.B,ge.Dv),[gx(k);gx(j)],act_z,[]);
                    end
                end
                pg (pg < max(pg(:))*0.01*probcutoff) = max(pg(:))*0.01*probcutoff;         
                imagesc(gx(cutoff:end),gx(cutoff:end),pg(cutoff:end,cutoff:end))   
                colormap jet
                hold on;
                plot(act_g(2),act_g(1),'-gx','MarkerSize',20,'LineWidth',3);

                subplot(1,4,2);
                imagesc(gx(cutoff:end),gx(cutoff:end),gg1(cutoff:end,cutoff:end))
                colormap jet
                maxval = max([ max(max(gg1(cutoff:end,cutoff:end))) abs( min(min(gg1(cutoff:end,cutoff:end))) ) ]);
                colorscale = 1e-5;
                caxis(colorscale*[-maxval maxval]);
                title('1')
                hold on;
                plot(act_g(2),act_g(1),'-gx','MarkerSize',20,'LineWidth',3);

                subplot(1,4,3);
                imagesc(gx(cutoff:end),gx(cutoff:end),gg2(cutoff:end,cutoff:end))
                colormap jet
                maxval = max([ max(max(gg2(cutoff:end,cutoff:end))) abs( min(min(gg2(cutoff:end,cutoff:end))) ) ]);                
                caxis(colorscale*[-maxval maxval]);
                title('2')
                hold on;
                plot(act_g(2),act_g(1),'-gx','MarkerSize',20,'LineWidth',3);      
                
                subplot(1,4,4);                
                zx = linspace(0,10,30);
                pz = zeros(length(zx),1);
                for j=1:length(zx)
                    pz(j) = gestaltFullLogPosterior(ge,act_x,act_v,act_g,zx(j),[]);
                end
                plot(zx,pz)                                

                pause
            end

            prev_g = act_g;
            prev_v = act_v;
            prev_z = act_z;

            max_shift_g = 0.1;
            act_v = act_v + act_lr_v * grad_V';
            act_g = max(act_g + min(act_lr_g * grad_G,max_shift_g),minval*ones(ge.k,1));
            act_z = max(act_z + act_lr_z * grad_Z,minval);            

            delta_g = max((prev_g-act_g).^2);
            delta_v = max((prev_v-act_v).^2);
            delta_z = max((prev_z-act_z).^2);
            maxdelta = max([delta_g delta_v delta_z]);
            if maxdelta < delta_threshold || counter == 5000
                convergence = true;
                fprintf('Convergence achieved in %d steps.\n',counter-1);
            elseif rem(counter-1,1000) == 0
                act_lr_g = act_lr_g * 0.1;
                act_lr_v = act_lr_v * 0.1;
                act_lr_z = act_lr_z * 0.1;
                delta_threshold = delta_threshold * 0.1;
            end

%             actnorm = reshape(act_v,1,ge.Dv);
%             actnorm = actnorm/norm(actnorm);
%             angle = actnorm * norm_v';
            
            if saveCourses
                actvcourse = [actvcourse act_v'];
                actgcourse = [actgcourse act_g];
                actzcourse = [actzcourse act_z];
                lp = gestaltFullLogPosterior(ge,act_x,act_v,act_g,act_z,[]);
                actloglike = [actloglike lp];
            end

        end
        
        if saveCourses
            vcourse{end+1} = actvcourse;
            gcourse{end+1} = actgcourse;
            zcourse{end+1} = actzcourse;
            loglike{end+1} = actloglike;
        end
        V(i,:) = act_v;
        G(i,:) = act_g';
        Z(i,1) = act_z;
    end
    
    if saveCourses
        V = vcourse;
        G = gcourse;
        Z = zcourse;
    end
end