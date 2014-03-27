function [prg,prv,pr] = gestaltPlotPosterior(ge,x)
    clf;
    if ge.Dv == 1
        prg = gestaltPlot1D(ge,x);       
        return
    elseif ge.Dv == 2
        return
    elseif ge.Dv == 3
        red = ge;
    else
        red = gestaltReduce(ge,2,3);
    end
    
    gmax = 1;
    gstep = 0.1;
    gx = 0.01:gstep:gmax-0.01;
    vx = -1:0.25:1;
    vy = -1:0.25:1;
    vz = -1:0.25:1;

    pr = zeros(size(gx,2),size(vx,2),size(vy,2),size(vz,2));
    %fprintf('Contour %d/', size(gx,2)*size(vx,2)*size(vy,2)*size(vz,2));
    for i=1:size(gx,2)
                            g = [gx(1,i);1-gx(1,i)];
         for k=1:size(vx,2)
             for l=1:size(vy,2)
                 for o=1:size(vz,2)
     %               printCounter((i-1)*size(vx,2)*size(vy,2)*size(vz,2) + (k-1)*size(vy,2)*size(vz,2) + (l-1)*size(vz,2) + o);                    
                    pr(i,k,l,o) = gestaltUnnormLogPost(g,[vx(1,k);vy(1,l);vz(1,o)],x,red);
                 end
             end
         end
    end
    fprintf('\n');
    
    
    prg = sum(sum(sum(pr,4),3),2);    
    prv = reshape(sum(sum(pr,1),2),size(vx,2),size(vy,2));

    subplot(1,2,1);
    plot(gx,prg);
    hold on;
    subplot(1,2,2);
    contour(vx,vy,prv);
    hold on;
end