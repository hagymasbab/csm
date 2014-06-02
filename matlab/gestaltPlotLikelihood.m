function ll = gestaltPlotLikelihood(ge,L,gridnum,complete,ll)
    % only works for 2-dimensional, 1-component models
    x = linspace(0,1,gridnum);
    y = linspace(0,1,gridnum);
    z = linspace(0,1,gridnum);
    tics = linspace(0,1,4);
    
    true_x = ge.cc{1}(1,1);
    true_y = ge.cc{1}(2,2);
    true_z = ge.cc{1}(1,2);
    [~,trueindex_x] = min(abs(x - true_x));
    [~,trueindex_y] = min(abs(y - true_y));
    [~,trueindex_z] = min(abs(z - true_z));
    
    if isscalar(L) && complete
        sdim = ge.k+(ge.Dv*ge.B);
        samples = zeros(ge.N,L,sdim);
        for n=1:ge.N
            [samples(n,:,:),~] = gestaltGibbs(ge,n,L,'verbose',0);
        end
    else
        samples = L;
    end
    
    if ll==0
        ll = zeros(gridnum,gridnum,gridnum);
        fprintf('%d/',gridnum^3);
        for i=1:gridnum
            for j=1:gridnum
                for k=1:gridnum
                    printCounter((i-1)*gridnum^2 + (j-1)*gridnum + k);
                    ge.cc{1} = [x(i) z(k); z(k) y(j)];
                    [~,err] = chol(ge.cc{1});
                    if err == 0 && rcond(ge.cc{1}) > 1e-16       
                        if complete
                            ch = cellchol(ge.cc);
                            ll(i,j,k) = gestaltCompleteDataLogLikelihood(ge,samples,ch);
                        else
                            ll(i,j,k) = gestaltLogLikelihood(ge,L,0,[]);
                        end
                    else
                        ll(i,j,k) = -Inf;
                    end
                end
            end
        end
        fprintf('\n');
        ll = shiftLogData(ll);
    end                        
    
    clf;
    subplot(4,3,1);
    viewImage(sum(ll,3),'usemax',true);
    xlabel('C(2,2)');
    ylabel('C(1,1)');

    subplot(4,3,2);
    viewImage(squeeze(sum(ll,2)),'usemax',true);
    xlabel('C(1,2)');
    ylabel('C(1,1)');
    
    subplot(4,3,3);
    viewImage(squeeze(sum(ll,1)),'usemax',true);
    xlabel('C(1,2)');
    ylabel('C(2,2)');      
    
    subplot(4,3,4);
    plot(x,sum(sum(ll,3),2));
    xlabel('C(1,1)');
    ylim = get(gca,'YLim');
    hold on
    plot([true_x true_x], [0 ylim(2)],'r')
    
    subplot(4,3,5);
    plot(y,sum(sum(ll,3),1));
    xlabel('C(2,2)');
    ylim = get(gca,'YLim');
    hold on
    plot([ge.cc{1}(2,2) ge.cc{1}(2,2)], [0 ylim(2)],'r')
    
    subplot(4,3,6);
    plot(z,reshape(sum(sum(ll,1),2),1,gridnum));
    xlabel('C(1,2)');
    ylim = get(gca,'YLim');
    hold on
    plot([ge.cc{1}(1,2) ge.cc{1}(1,2)], [0 ylim(2)],'r')
    
    subplot(4,3,7);
    viewImage(squeeze(ll(:,:,trueindex_z)),'usemax',true);
    xlabel('C(2,2)');
    ylabel('C(1,1)');
    title(sprintf('C(1,2) = %.2f',z(trueindex_z)));
    
    subplot(4,3,8);
    viewImage(squeeze(ll(:,trueindex_y,:)),'usemax',true);
    xlabel('C(1,2)');
    ylabel('C(1,1)');
    title(sprintf('C(2,2) = %.2f',y(trueindex_y)));
    
    subplot(4,3,9);
    viewImage(squeeze(ll(trueindex_x,:,:)),'usemax',true);
    xlabel('C(1,2)');
    ylabel('C(2,2)');
    title(sprintf('C(1,1) = %.2f',x(trueindex_x)));
    
    subplot(4,3,10);
    plot(x,squeeze(ll(:,trueindex_y,trueindex_z))');
    xlabel('C(1,1)');
    ylim = get(gca,'YLim');
    hold on
    plot([true_x true_x], [0 ylim(2)],'r')
    title(sprintf('C(2,2)=%.2f C(1,2)=%.2f',y(trueindex_y),z(trueindex_z)));
    
    subplot(4,3,11);
    plot(y,squeeze(ll(trueindex_x,:,trueindex_z))');
    xlabel('C(2,2)');
    ylim = get(gca,'YLim');
    hold on
    plot([true_y true_y], [0 ylim(2)],'r')
    title(sprintf('C(1,1)=%.2f C(1,2)=%.2f',x(trueindex_x),z(trueindex_z)));
    
    subplot(4,3,12);
    plot(z,squeeze(ll(trueindex_x,trueindex_y,:))');
    xlabel('C(1,2)');
    ylim = get(gca,'YLim');
    hold on
    plot([true_z true_z], [0 ylim(2)],'r')
    title(sprintf('C(1,1)=%.2f C(2,2)=%.2f',x(trueindex_x),y(trueindex_y)));

end
    
function md = mindiff(A)
    md = realmax;
    for i=1:numel(A)
        for j=i+1:numel(A)
            if abs(i-j)> 0 && abs(i-j) < md
                md = abs(i-j);
            end
        end
    end
end

function sd = shiftLogData(d)
    md = mindiff(d);
    noinf = d(d~=-Inf);
    minel = min(noinf(:));
    sd = d;
    sd(d==-Inf) = minel-md;
    sd = sd - (minel-md);
end