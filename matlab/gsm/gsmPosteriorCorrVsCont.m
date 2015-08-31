function gsmPosteriorCorrVsCont(A,C,sigma_x,zs)
    nZ = length(zs);
    iC = stableInverse(C);
    sATA = (1/sigma_x^2) * (A' * A);
    ymax = 0;
    for i=1:nZ
        actCov = stableInverse(iC + zs(i)^2 * sATA);
        actCorr = upperTriangleValues(corrcov(actCov));
        subplot(1,nZ,i);
        hist(actCorr,linspace(-1,1,100));
        yl = ylim();
        if yl(2) > ymax
            ymax = yl(2);
        end
    end
    for i=1:nZ
        subplot(1,nZ,i);
        xlim([-1 1]);
        ylim([0 ymax*1.05]);
    end
end