function compareCorrelations(dim,filters,t1,t2)

    filterfile = sprintf('filters_%s_%d.mat',filters,dim);

    [cm1,U1] = linearFilterCorrelation(1000,sprintf('patches_text%d_%d.mat',t1,dim),filterfile);
    [cm2,U2] = linearFilterCorrelation(1000,sprintf('patches_text%d_%d.mat',t2,dim),filterfile);
    diff12 = cm1 - cm2;
    [~,idx] =max(diff12(:));
    [a,b] = ind2sub(size(cm1),idx);

    load(filterfile);
    subplot(3,2,1);
    viewImage(A(:,a),'useMax',true);
    title('Filter 1','FontSize',16);
    subplot(3,2,2);
    viewImage(A(:,b),'useMax',true);
    title('Filter 2','FontSize',16);
    
    subplot(3,2,3);
    scatter(U1(:,a),U1(:,b),'k');
    yl1 = ylim();
    xl1 = xlim();
    title('Responses for texture type 1','FontSize',16);
    subplot(3,2,4);
    scatter(U2(:,a),U2(:,b),'k');
    yl2 = ylim();
    xl2 = xlim();
    title('Responses for texture type 2','FontSize',16);
    subplot(3,2,3);
    xlim([ min([xl1(1) xl2(1)]) max([xl1(2) xl2(2)]) ]);
    ylim([ min([yl1(1) yl2(1)]) max([yl1(2) yl2(2)]) ]);
    hold on
    correlationPlot(U1(:,a),U1(:,b),cm1(a,b));
    hold off
    subplot(3,2,4);
    xlim([ min([xl1(1) xl2(1)]) max([xl1(2) xl2(2)]) ]);
    ylim([ min([yl1(1) yl2(1)]) max([yl1(2) yl2(2)]) ]);
    hold on
    correlationPlot(U2(:,a),U2(:,b),cm2(a,b));
    hold off
    
    subplot(3,2,5);
    imgsamp1 = imread(sprintf('text%d_sample.jpg',t1));
    viewImage(imgsamp1,'useMax',true);
    title('Texture type 1','FontSize',16);
    subplot(3,2,6);
    imgsamp2 = imread(sprintf('text%d_sample.jpg',t2));
    viewImage(imgsamp2,'useMax',true);
    title('Texture type 2','FontSize',16);
end