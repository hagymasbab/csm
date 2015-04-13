function plotMAPresults(v,g,z,loglike,ge)

    close all
    colnum = 5;    

    for ni = 1:length(v)
        figure

        map_v = v{ni}(:,end);
        init_v = v{ni}(:,1);
        map_z = z{ni}(end,1);
        stepnum = length(z{ni});        
        v_true = reshape(ge.V(ni,1,:),ge.Dv,1);
        length_true = norm(v_true);
        norm_true = v_true / length_true; 
        
        angle = zeros(stepnum,1);
        normdist = zeros(stepnum,1);
        distance = zeros(stepnum,1);
        for i=1:stepnum
            v_act = v{ni}(:,i);
            length_act = norm(v_act);
            norm_act = v_act/length_act;
            
            angle(i,1) = norm_act' * norm_true;            
            normdist(i,1) = (length_act-length_true)^2;
            distance(i,1) = sum((v_act - v_true).^2);
        end
        
        
        subplot(2,colnum,1);
        viewImage(map_v);
        title('MAP V');

        subplot(2,colnum,2);
        viewImage(ge.V(ni,1,:));
        title('True V');

        subplot(2,colnum,colnum+1);
        viewImage(map_z * ge.A * map_v + ge.obsVar*randn(ge.Dv,1));
        title('X from MAP V and Z');

        subplot(2,colnum,colnum+2);
        viewImage(ge.X(ni,1,:));
        title('True X');

        subplot(2,colnum,3);
        %viewImage(ge.A'*reshape(ge.X(ni,1,:),1,ge.Dv)');
        viewImage(init_v);
        title('V init');

        subplot(2,colnum,colnum+3);
        ymax = max([1 max(distance) max(normdist) max(angle)]);
        plot(angle*ymax,'b','LineWidth',2);
        hold on;
        plot(normdist,'g','LineWidth',2);
        plot(distance,'r','LineWidth',2);
        xlim([1 stepnum]);
        ylim([0 ymax + 0.1]);
        title('V_{MAP} \cdot V_{True}');

        subplot(2,colnum,4);
        plot(g{ni}(1,:));
        max1 = max(g{ni}(1,:));
        hold on
        plot([0;stepnum],[ge.G(ni,1) ge.G(ni,1)],'r-');
        title('G 1');
        xlim([1 stepnum]);

        subplot(2,colnum,colnum+4);
        plot(g{ni}(2,:));
        hold on
        plot([0;stepnum],[ge.G(ni,2) ge.G(ni,2)],'r-');
        title('G 2');
        xlim([1 stepnum]);
        max2 = max(g{ni}(2,:));
        ymax = max([max1 max2 ge.G(ni,1) ge.G(ni,2)]) + 0.1;
        ylim([0 ymax]);
        subplot(2,colnum,4);
        ylim([0 ymax]);

        subplot(2,colnum,5);
        plot(z{ni});
        hold on
        plot([0;stepnum],[ge.Z(ni,1) ge.Z(ni,1)],'r-');
        title('Z');
        xlim([1 stepnum]);

        subplot(2,colnum,colnum+5);
        plot(loglike{ni});    
        title('Log-likelihood');
        xlim([1 stepnum]);
    end
end