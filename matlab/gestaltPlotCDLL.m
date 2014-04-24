function ll=gestaltPlotCDLL(ge,samples,gridnum,negative,ll)
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
                        ch = cellchol(ge.cc);
                        ll(i,j,k) = exp(gestaltCompleteDataLogLikelihood(ge,samples,ch));
                    else
                        ll(i,j,k) = 0;
                    end
                end
            end
        end
        fprintf('\n');
    end        
    
    subplot(2,3,1);
    imagesc(sum(ll,3));
    xlabel('C(2,2)');
    ylabel('C(1,1)');
    set(gca,'XTick',tics);
    set(gca,'YTick',tics);
    subplot(2,3,2);
    imagesc(squeeze(sum(ll,2)));
    xlabel('C(1,2)');
    ylabel('C(1,1)');
    subplot(2,3,3);
    imagesc(squeeze(sum(ll,1)));
    xlabel('C(1,2)');
    ylabel('C(2,2)');
    subplot(2,3,4);
    plot(x,sum(sum(ll,3),2));
    xlabel('C(1,1)');
    ylim = get(gca,'YLim');
    line([ge.cc{1}(1,1) ge.cc{1}(1,1)], [0 ylim(2)])
    subplot(2,3,5);
    plot(y,sum(sum(ll,3),1));
    xlabel('C(2,2)');
    ylim = get(gca,'YLim');
    line([ge.cc{1}(2,2) ge.cc{1}(2,2)], [0 ylim(2)])
    subplot(2,3,6);
    plot(z,reshape(sum(sum(ll,1),2),1,gridnum));
    xlabel('C(1,2)');
    ylim = get(gca,'YLim');
    line([ge.cc{1}(1,2) ge.cc{1}(1,2)], [0 ylim(2)])
    