function compareCorrelations(nPatches,dim,filters,t1,t2,gsm_res,predef_pair,loadStuff)

    filterfile = sprintf('filters_%s_%d.mat',filters,dim);

    if loadStuff
        load('save_compcorr.mat');
    else
        [cm1,U1] = linearFilterCorrelation(nPatches,sprintf('patches_text%d_%d.mat',t1,dim),filterfile,gsm_res);
        [cm2,U2] = linearFilterCorrelation(nPatches,sprintf('patches_text%d_%d.mat',t2,dim),filterfile,gsm_res);
        if gsm_res > 0
            save('bin/save_compcorr.mat','cm1','cm2','U1','U2');
        end
    end
    
    if isempty(predef_pair)
        diff12 = cm1 - cm2;
        [~,idx] =max(diff12(:));
        [a,b] = ind2sub(size(cm1),idx);
    else
        a = predef_pair(1);
        b = predef_pair(2);
    end

    load(filterfile);
    subplot(3,2,1);
    viewImage(A(:,a),'useMax',true);
    title(sprintf('Filter %d',a),'FontSize',16);
    subplot(3,2,2);
    viewImage(A(:,b),'useMax',true);
    title(sprintf('Filter %d',b),'FontSize',16);
    
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