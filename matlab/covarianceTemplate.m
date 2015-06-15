function template = covarianceTemplate(filters,aspects,thresholds)   
    close all;    
    dv = size(filters,2);
    template = false(dv);
    plot_num = 0;
    plot_count = 0;
    to_plot = cell(plot_num,2);
    
    orients = zeros(dv,1);    
    for i = 1:dv
        orients(i) = orientation(filters(:,i));
    end
    
    for a=1:length(aspects)      
        subtemplate = false(dv);
        for i = 1:dv
            for j = i+1:dv
                subtemplate(i,j) = feval(aspects{a},i,j,thresholds{a},filters,orients);
                if subtemplate(i,j) && plot_count < plot_num    
                    plot_count = plot_count + 1;
                    to_plot{plot_count,1} = filters(:,i);
                    to_plot{plot_count,2} = filters(:,j);
                end
            end
        end
        template = template | subtemplate;
    end
    template = template + logical(eye(dv));
    if plot_num > 0
        viewImageSet(to_plot);
        figure();
        viewImage(template);
    end
end

function passed = overlap(i1,i2,threshold,filters,orients) 
    f1 = filters(:,i1);
    f2 = filters(:,i2);
    passed = abs(f1' * f2) > threshold;        
end

function passed = parallell(i1,i2,threshold,filters,orients)
    o1 = orients(i1);
    o2 = orients(i2);
    passed = abs(o1-o2) < threshold || abs(o1-180 - o2) < threshold || abs(o1 - o2-180) < threshold;
end

function theta = orientation(f)
    dv = length(f);
    imdim = sqrt(dv);
    [~,maxes] = max(f);
    maxX = ceil(maxes / imdim);
    maxY = rem(maxes,imdim);
    
    th_range = 1:180;
    maxprod = 0;
    theta = -1;
    prods = zeros(length(th_range));
    for i=1:length(th_range)
        a = tand(th_range(i));
        b = maxY - maxX * a;
        pos = ones(imdim);
        for j = 1:imdim
            for k = 1:imdim
                if j < floor(k*a+b)
                    pos(j,k) = -1;
                end
            end
        end        
        neg = -pos;
        actprod = max(f' * pos(:),f' * neg(:));
        prods(i) = actprod;
        if  actprod > maxprod
            maxprod = actprod;
            theta = th_range(i);
        end
        %subplot(1,3,1);viewImage(pos);subplot(1,3,2);viewImage(f,'useMax',true);subplot(1,3,3);plot(prods);pause;
    end
end