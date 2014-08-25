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
    gx = 0.01:0.01:0.99;
    gy = 0.01:0.01:0.99;
    
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
        for g=1:size(gx,2)
            if strcmp(prior,'dirichlet')
                lp(1,g) = gestaltLogPostG(gx(g),actV,ge,prior,precision);
            else
                for gg=1:size(gy,2)
                    lp2dim(g,gg) = gestaltLogPostG([gx(g);gy(gg)],actV,ge,prior,precision);
                end
            end                
        end

        if ~strcmp(prior,'dirichlet')
            lp1 = sum(lp2dim,2);
            lp2 = sum(lp2dim,1);

        else
            lp1 = lp;
            lp2 = ones(1,size(gx,2)) - lp1; 
        end
        
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