function plotTuningCurves()
    load('gsm_tc_newlearn1.mat')
    tc_early = tuningCurves;
    load('gsm_tc_newlearn19.mat')
    tc_late = tuningCurves;
    
    %to_plot = [1 2 3 4 5];    
    %to_plot = chooseKfromN(10,size(tc_early,1));
    to_plot = [5 111 28 144 241];
    nPlot = length(to_plot);
    ymax = 0;
    for i=1:nPlot
        subplot(nPlot,2,(i-1)*2+1);
        plot(linspace(1,180,18),tc_early(to_plot(i),:));        
        hold on;
        plot(linspace(1,180,18),tc_early(to_plot(i),:),'LineWidth',3);
        hold off;
        xlim([0,180]);
        yl = ylim();
        if yl(2) > ymax
            ymax = yl(2);
        end
        set(gca,'FontSize',16,'Ytick',[]);
        if i<nPlot
            set(gca,'Xtick',[]);
        end
            
        subplot(nPlot,2,i*2);
        plot(linspace(1,180,18),tc_late(to_plot(i),:),'LineWidth',3);
        xlim([0,180]);
        yl = ylim();
        if yl(2) > ymax
            ymax = yl(2);
        end
        set(gca,'FontSize',16,'Ytick',[]);
        if i<nPlot
            set(gca,'Xtick',[]);
        end
    end
    ymax = 6;
    for i=1:nPlot
        subplot(nPlot,2,(i-1)*2+1);
        ylim([0 ymax]);
        subplot(nPlot,2,i*2);
        ylim([0 ymax]);
    end
end