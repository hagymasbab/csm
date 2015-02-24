function bcc = parametrisationDifference(ge,cc,loaddata)
    if loaddata
        load('bin/save_paramdiff.mat');
    else
        ccm = cell(1,length(cc));
        bcm = cell(1,length(cc));
        bcc = cell(1,length(cc));
        for i=1:length(cc)
            b = sign(cc{i}(1,:))' .* sqrt(diag(cc{i}));
            bcc{i} = b*b';
            ccm{i} = corrcov(cc{i});
            bcm{i} = corrcov(bcc{i});
        end
    %     viewImageSet(bcc)

        [~,seeds,~,cst] = gestaltGReceptiveFields(ge,cc,10000,false);
        [~,bseeds,~,bst] = gestaltGReceptiveFields(ge,bcc,10000,false);
        save('bin/save_paramdiff.mat','bcc','seeds','bseeds','cst','bst','cc','ccm','bcm');
    end
    hornum = 6;
    subplot = @(m,n,p) subtightplot (m, n, p, [0.025 0.001], [0 0.025], [0 0.01]);
    load cmp_diff
    for i=1:length(cc)      
        printCounter(i,'maxVal',length(cc),'stringVal','Component');
        diffcm = ccm{i} - bcm{i};
        diffst = cst{i} - bst{i};
        
        um = true;
        subplot(length(cc),hornum,(i-1)*hornum + 1);
        viewImage(ccm{i},'useMax',um);
        if i==1;title('Corr.mat.');end;
        ylabel(sprintf('Component %d',i));
        freezeColors
        subplot(length(cc),hornum,(i-1)*hornum + 2);
        viewImage(bcm{i},'useMax',um);
        if i==1;title('Corr.mat. KL');end;
        freezeColors
        subplot(length(cc),hornum,(i-1)*hornum + 3);
        viewImage(diffcm,'useMax',um);
        if i==1;title('Corr.mat. Diff');end;
        colormap(cmp_diff);
        freezeColors
        
        subplot(length(cc),hornum,(i-1)*hornum + 4);
        viewImage(cst{i},'useMax',um);
        if i==1;title('Pixel');end;
        freezeColors
        subplot(length(cc),hornum,(i-1)*hornum + 5);
        viewImage(bst{i},'useMax',um);
        if i==1;title('Pixel KL');end;
        freezeColors
        subplot(length(cc),hornum,(i-1)*hornum + 6);
        viewImage(diffst,'useMax',um);
        if i==1;title('Pixel Diff');end;
        colormap(cmp_diff);
        freezeColors
    end
end
        
        