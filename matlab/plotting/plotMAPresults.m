close all
colnum = 5;

for ni = 1:2
    figure
    
    subplot(2,colnum,1);
    viewImage(v(ni,:));
    title('MAP V');

    subplot(2,colnum,2);
    viewImage(ge.V(ni,1,:));
    title('True V');

    subplot(2,colnum,colnum+1);
    viewImage(z(ni,1) * ge.A * v(ni,:)' + ge.obsVar*randn(ge.Dv,1));
    title('X from MAP V and Z');

    subplot(2,colnum,colnum+2);
    viewImage(ge.X(ni,1,:));
    title('True X');

    subplot(2,colnum,3);
    viewImage(ge.A'*reshape(ge.X(ni,1,:),1,ge.Dv)');
    title('V init');

    subplot(2,colnum,colnum+3);
    plot(delta{ni}');
    xlim([1 length(delta{ni})]);
    title('V_{MAP} \cdot V_{True}');

    subplot(2,colnum,4);
    plot(gcourse{ni}(1,:));
    hold on
    plot([0;length(delta{ni})],[ge.G(ni,1) ge.G(ni,1)],'r-');
    title('G 1');
    xlim([1 length(delta{ni})]);

    subplot(2,colnum,colnum+4);
    plot(gcourse{ni}(2,:));
    hold on
    plot([0;length(delta{ni})],[ge.G(ni,2) ge.G(ni,2)],'r-');
    title('G 2');
    xlim([1 length(delta{ni})]);

    subplot(2,colnum,5);
    plot(zcourse{ni});
    hold on
    plot([0;length(delta{ni})],[ge.Z(ni,1) ge.Z(ni,1)],'r-');
    title('Z');
    xlim([1 length(delta{ni})]);
    
    subplot(2,colnum,colnum+5);
    plot(loglike{ni});    
    title('Log-likelihood');
    xlim([1 length(delta{ni})]);
end