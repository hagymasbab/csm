function posteriors = gestaltPlotCondPostG(ge,V,prior,stepsize)
    %clf;
    precision = false;
    nInput = 1;
    if ndims(V)>2
        if size(V,1) == 1
            V = reshape(V,ge.B,ge.Dv);
        else
            nInput = size(V,1);
        end
    end
        
    % k=2
    %stepsize = 0.1;
    gmax = 5;
    gx = 0.01:stepsize:gmax;
    gy = 0.01:stepsize:gmax;
    gz = 0.01:stepsize:gmax;
    
    legends = cell(1,nInput);
    maxindices = zeros(1,nInput);
    colors = {'b','g','r'};
    handles = zeros(1,nInput);
    for i = 1:nInput
        if nInput > 1
            actV = reshape(V(i,:,:),ge.B,ge.Dv);
        else
            actV = V;
        end
        lp = zeros(1,size(gx,2));
        lp2dim = zeros(size(gx,2),size(gy,2));
        lp3dim = zeros(size(gx,2),size(gy,2),size(gz,2));
        pr2dim = zeros(size(gx,2),size(gy,2));
        pr3dim = zeros(size(gx,2),size(gy,2),size(gz,2));
        for g=1:size(gx,2)
            printCounter(g,'stringVal','first dim','maxVal',size(gx,2));
            if strcmp(prior,'dirichlet')
                lp(1,g) = gestaltLogPostG(gx(g),actV,ge,prior,precision);
            else                
                for gg=1:size(gy,2)
                    if ge.k == 2
                        lp2dim(g,gg) = exp(gestaltLogPostG([gx(g);gy(gg)],actV,ge,prior,precision));
                        pr2dim(g,gg) = exp(gestaltLogPriorG([gx(g);gy(gg)],ge));
                    elseif ge.k >= 3
                        for ggg = 1:size(gz,2)
                            if ge.k == 4
                                actlp = 0;
                                actpr = 0;
                                for g4 = 1:size(gz,2)
                                    actlp = actlp + gestaltLogPostG([gx(g);gy(gg);gz(g4);gz(ggg)],actV,ge,prior,precision);
                                    actpr = actpr + gestaltLogPriorG([gx(g);gy(gg);gz(g4);gz(ggg)],ge);
                                end
                                lp3dim(g,gg,ggg) = exp(actlp);
                                pr3dim(g,gg,ggg) = exp(actpr);
                            else
                                lp3dim(g,gg,ggg) = exp(gestaltLogPostG([gx(g);gy(gg);gz(ggg)],actV,ge,prior,precision));
                                pr3dim(g,gg,ggg) = exp(gestaltLogPriorG([gx(g);gy(gg);gz(ggg)],ge));
                            end
                        end
                    end                
                end
            end                
        end

        if strcmp(prior,'gamma')
            if ge.k == 2
                lp1 = sum(lp2dim,2);
                lp2 = sum(lp2dim,1);
                pr1 = sum(pr2dim,2);
                pr2 = sum(pr2dim,1);
                lp12 = lp2dim;
                pr12 = pr2dim;
                legends = {'comp1post','comp0post'};
            elseif ge.k >= 3
                lp12 = squeeze(sum(lp3dim,3));
                lp13 = squeeze(sum(lp3dim,2));
                %lp23 = squeeze(sum(lp3dim,1));
                lp1 = sum(lp12,2);
                lp2 = sum(lp12,1);
                lp3 = sum(lp13,1);
                pr12 = squeeze(sum(pr3dim,3));
                pr1 = squeeze(pr3dim(:,1,1));
                pr2 = squeeze(pr3dim(1,:,1));
                pr3 = squeeze(pr3dim(1,1,:));
                legends = {'comp1post','comp2post','comp0post'};
            end
            
            lp1 = lp1 / norm(lp1);
            lp2 = lp2 / norm(lp2);            
            pr1 = pr1 / norm(pr1);
            pr2 = pr2 / norm(pr2);
            posteriors = [lp1(:)';lp2(:)'];                
            priors = [pr1(:)';pr2(:)'];
            
            if ge.k >= 3
                lp3 = lp3 / norm(lp3);                
                pr3 = pr3 / norm(pr3);
                posteriors = [posteriors;lp3(:)'];                
                priors = [priors;pr3(:)'];                
            end
                 
            subplot(1,2,1)
            posteriors = shiftSignals(posteriors);
            plot(gx,posteriors);
            hold on;                
            priors = shiftSignals(priors);                
            plot(gx,priors,'--');                            
            legend(legends);
            hold off;
                        
            subplot(1,2,2)
            
            hold on;
            contour(gx,gy,lp12);
            contour(gx,gy,pr12);
                        
            return;
            % TODO this is disgusting
            
        else
            lp1 = lp;
            lp2 = ones(1,size(gx,2)) - lp1; 
        end
        
        %if ge.k == 3
            
            
        [~,maxindices(i)] = max(lp1);
        handles(i) = plot(gx,lp1,colors{i});
        hold on;
        xval = gx(maxindices(i));
        plot([xval xval], ylim, colors{i});
        if nInput == 1
            plot(gy,lp2,'r');
        else
            legends{i} = sprintf('%d',i);
%             legends{nInput + i} = sprintf('max%d',i);
        end
    end
    
%     for i = 1:nInput
%         xval = gx(maxindices(i));
%         plot([xval xval], ylim);
%         hold all;
%     end
%   
    if nInput > 1
        legend(handles,legends);
    end
    hold off;
end