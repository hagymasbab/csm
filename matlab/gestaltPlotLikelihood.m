function ll = gestaltPlotLikelihood(ge,L,gridnum,negative,ll)
    % only works for 2-dimensional, 1-component models
    x = linspace(-1*double(negative),1,gridnum);
    y = linspace(-1*double(negative),1,gridnum);
    z = linspace(-1*double(negative),1,gridnum);
    tics = linspace(-1*double(negative),1,4);
    
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
                        [co,ex] = gestaltLikelihood(ge,L);
                        ll(i,j,k) = log10(co) + ex;
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
    subplot(2,3,1);
    viewImage(sum(ll,3),'usemax',true);
    xlabel('C(2,2)');
    ylabel('C(1,1)');
    set(gca,'XTick',tics);
    set(gca,'YTick',tics);
    subplot(2,3,2);
    viewImage(squeeze(sum(ll,2)),'usemax',true);
    %imagesc(squeeze(sum(ll,2)));
    xlabel('C(1,2)');
    ylabel('C(1,1)');
    subplot(2,3,3);
    viewImage(squeeze(sum(ll,1)),'usemax',true);
    %imagesc(squeeze(sum(ll,1)));
    xlabel('C(1,2)');
    ylabel('C(2,2)');      
    
    subplot(2,3,4);
    plot(x,sum(sum(ll,3),2));
    xlabel('C(1,1)');
    ylim = get(gca,'YLim');
    hold on
    plot([ge.cc{1}(1,1) ge.cc{1}(1,1)], [0 ylim(2)],'r')
    
    subplot(2,3,5);
    plot(y,sum(sum(ll,3),1));
    xlabel('C(2,2)');
    ylim = get(gca,'YLim');
    hold on
    plot([ge.cc{1}(2,2) ge.cc{1}(2,2)], [0 ylim(2)],'r')
    
    subplot(2,3,6);
    plot(z,reshape(sum(sum(ll,1),2),1,gridnum));
    xlabel('C(1,2)');
    ylim = get(gca,'YLim');
    hold on
    plot([ge.cc{1}(1,2) ge.cc{1}(1,2)], [0 ylim(2)],'r')

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