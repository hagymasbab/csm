function findSameVar(cc,loaddata)
    Dx = size(cc{1},1);
    k = length(cc);
    numel = 2;
     
    if loaddata
        load('bin/save_samevar.mat');
    else
        pairvardiffs = zeros(Dx,Dx,k,k);
        for i = 1:Dx
            for j=i+1:Dx
                for ik = 1:k
                    for jk = ik+1:k
                        pairvardiffs(i,j,ik,jk) = abs(cc{ik}(i,i) - cc{jk}(i,i)) + abs(cc{ik}(j,j) - cc{jk}(j,j));
                    end
                end
            end
        end

        cm = cell(1,length(cc));
        for i=1:length(cc)        
            cm{i} = corrcov(cc{i});
        end

        pairvardiffs(pairvardiffs==0) = Inf;
        [~,sortIndex] = sort(pairvardiffs(:),'ascend');
        minIndex = sortIndex(1:numel);    
        save('bin/save_samevar.mat','pairvardiffs','cm','minIndex')
    end
        
    %hist(pairvardiffs(:),200)
   
    
    horplot = 2;
    varmax = 0;
    corrmin = 0;
    corrmax = 0;
    load('cmp_graybars.mat');
    colormap(cmp_graybars);
    for i=1:numel
        [cell1,cell2,comp1,comp2] = ind2sub(size(pairvardiffs),minIndex(i));        
        fprintf('Cells %d and %d\n',cell1,cell2);
        fprintf('\tin component %d: variance1 %.2f variance2 %.2f correlation %.2f\n',comp1,cc{comp1}(cell1,cell1),cc{comp1}(cell2,cell2),cm{comp1}(cell1,cell2));
        fprintf('\tin component %d: variance1 %.2f variance2 %.2f correlation %.2f\n',comp2,cc{comp2}(cell1,cell1),cc{comp2}(cell2,cell2),cm{comp2}(cell1,cell2));
        
        subplot(numel,horplot,(i-1)*horplot+1);
        barh([cc{comp1}(cell1,cell1) cc{comp1}(cell2,cell2); cc{comp2}(cell1,cell1) cc{comp2}(cell2,cell2)]);
        ylim([0.5 2.5]);
        set(gca,'YtickLabel',{sprintf('%d',comp1),sprintf('%d',comp2)},'FontSize',16);
        if i<numel;set(gca,'Xtick',[]);end;
        xlim([0 max([cc{comp1}(cell1,cell1);cc{comp1}(cell2,cell2)])+0.02]);
        xlims = xlim();
        if xlims(2) > varmax; varmax = xlims(2);end;
        
        subplot(numel,horplot,(i-1)*horplot+2);
        barh([cm{comp1}(cell1,cell2);cm{comp2}(cell1,cell2)]);
        ylim([0.5 2.5]);        
        set(gca,'Ytick',[],'FontSize',16);
        if i<numel;set(gca,'Xtick',[]);end;
        xlim([min([cm{comp1}(cell1,cell2);cm{comp2}(cell1,cell2)])-0.02 max([cm{comp1}(cell1,cell2);cm{comp2}(cell1,cell2)])+0.02]);
        xlims = xlim();
        if xlims(2) > corrmax; corrmax = xlims(2);end;
        if xlims(1) < corrmin; corrmin = xlims(1);end;
    end
    
    for i=1:numel
        subplot(numel,horplot,(i-1)*horplot+1);
        xlim([0 varmax]);
        subplot(numel,horplot,(i-1)*horplot+2);
        xlim([corrmin corrmax]);        
    end                
end