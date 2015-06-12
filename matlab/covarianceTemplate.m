function template = covarianceTemplate(filters,aspects,thresholds)   
    close all;
    dv = size(filters,2);
    template = false(dv);
    plot_num = 10;
    plot_count = 0;
    to_plot = cell(plot_num,2);
    for a=1:length(aspects)      
        subtemplate = false(dv);
        for i = 1:dv
            for j = i+1:dv
                subtemplate(i,j) = feval(aspects{a},filters(:,i),filters(:,j),thresholds{a});
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
    viewImageSet(to_plot);
    figure();
    viewImage(template);
end

function passed = overlap(f1,f2,threshold) 
    passed = abs(f1' * f2) > threshold;        
end

function subtemplate = parallell(filters)
end