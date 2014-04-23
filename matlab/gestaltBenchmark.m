function diffs = gestaltBenchmark(ge,N,nRun,nSamples,maxStep,name)
    diffs = zeros(nRun,maxStep+1);
    fprintf('Run %d/',nRun);
    for r=1:nRun
        printCounter(r);
        ge = gestaltGenerate(ge,N,'verbose',false);
        diffs(r,:) = gestaltIEM(ge,ge.X,nSamples,maxStep,'shuffle','plot',0,'verbose',1);
        save(sprintf('diffs_%s.mat',name),'diffs');
    end
    fprintf('\n');
    %h=figure('visible','off');
    h=figure();
    errorbar(0:size(diffs,2)-1,mean(diffs,1),std(diffs,1));
    saveas(h,sprintf('convplot_%s.fig',name),'fig');
end