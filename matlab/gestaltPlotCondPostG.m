function lp = gestaltPlotCondPostG(ge,V,prior,precision)
    %clf;
    nInput = 1;
    if ndims(V)>2
        if size(V,1) == 1
            V = reshape(V,ge.B,ge.Dv);
        else
            nInput = size(V,1);
        end
    end
        
    % k=2
    gx = 0.01:0.05:0.99;
    gy = 0.01:0.05:0.99;
    gz = 0.01:0.05:0.99;
    
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
        pr3dim = zeros(size(gx,2),size(gy,2),size(gz,2));
        for g=1:size(gx,2)
            if strcmp(prior,'dirichlet')
                lp(1,g) = gestaltLogPostG(gx(g),actV,ge,prior,precision);
            else                
                for gg=1:size(gy,2)
                    if ge.k == 2
                        lp2dim(g,gg) = gestaltLogPostG([gx(g);gy(gg)],actV,ge,prior,precision);
                    elseif ge.k == 3
                        for ggg = 1:size(gz,2)
                             lp3dim(g,gg,ggg) = exp(gestaltLogPostG([gx(g);gy(gg);gz(ggg)],actV,ge,prior,precision));
                             pr3dim(g,gg,ggg) = exp(gestaltLogPriorG([gx(g);gy(gg);gz(ggg)],ge,prior));
                        end
                    end                
                end
            end                
        end

        if ~strcmp(prior,'dirichlet')
            if ge.k == 2
                lp1 = sum(lp2dim,2);
                lp2 = sum(lp2dim,1);
            elseif ge.k == 3
                lp12 = squeeze(sum(lp3dim,3));
                lp13 = squeeze(sum(lp3dim,2));
                %lp23 = squeeze(sum(lp3dim,1));
                lp1 = sum(lp12,2);
                lp2 = sum(lp12,1);
                lp3 = sum(lp13,2);
                posteriors = [lp1(:)';lp2(:)';lp3(:)'];
                
                step = (max(posteriors(:)) - min(posteriors(:))) / 100;
                
                posteriors = posteriors + repmat((0:2)*step,size(gx,2),1)';
                
                plot(gx,posteriors);
                hold on;
                
                pr12 = squeeze(sum(pr3dim,3));
                pr13 = squeeze(sum(pr3dim,2));                
                pr1 = sum(pr12,2);
                pr2 = sum(pr12,1);
                pr3 = sum(pr13,2);         
                priors = [pr1(:)';pr2(:)';pr3(:)'];
                
                priors = priors + repmat((0:2)*step,size(gx,2),1)';
                
                plot(gx,priors,'--');
                
                % TODO: nem jok itt a szinek!
                legend({'comp1post','comp2post','comp0post'})
                
                
                hold off;
                return;
                % TODO this is disgusting
            end
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