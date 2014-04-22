function gestaltBenchmark(ge,X,nRun,nSamples,maxStep,name)
    diffs = zeros(nRun,maxStep+1);
    for r=1:nRun
        diffs(r,:) = gestaltIEM(ge,X,nSamples,maxStep,'shuffle','plot',0);
        save(sprintf('diffs_%s.mat',name),'diffs');
    end
    %h=figure('visible','off');
    h=figure();
    errorbar(0:size(diffs,2)-1,mean(diffs,1),std(diffs,1));
    saveas(h,sprintf('convplot_%s.fig',name),'fig');
end