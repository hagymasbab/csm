function gestaltBowtie(nSamp,loadfile)
    close all;
    picNum = 9;
    delta_x = 0.1;
    plotnum = 12;
    Dx = 1024;
    firstPlotted = chooseKfromN(plotnum,Dx);
    fprintf('With contrast\n');
    plotSamples(picNum,delta_x,firstPlotted,true,loadfile,nSamp,plotnum,Dx);
    figure;
    fprintf('Without contrast\n');
    plotSamples(picNum,delta_x,firstPlotted,false,loadfile,nSamp,plotnum,Dx);
end

function plotSamples(picNum,delta_x,firstPlotted,contrast,loadfile,nSamp,plotnum,Dx)

    contstring = '';
    if contrast
        contstring = '_cont';
    end
    model_name = sprintf('jaguar%s',contstring);
    filename = sprintf('gestalt_%s.mat',model_name);
    if loadfile && exist(filename,'file') == 2
        load(filename);
        ge = gestalt;
    else
        ge = gestaltCreate(model_name,'Dx',Dx,'B',10,'obsVar',delta_x,'filters','gabor','N',picNum,'contrast',contrast);
    end
    s = zeros(picNum,nSamp,ge.k+ge.B*ge.Dx);
    zs = zeros(picNum,nSamp);
    ge.X = gestaltImageStimulus(picNum,ge.B,delta_x);
    fprintf('Sampling for image %d/',picNum);
    for p=1:picNum
        printCounter(p);
        [s(p,:,:),~,zs(p,:)] = gestaltGibbs(ge,p,nSamp,'contrast',contrast);
    end
    fprintf('\n');
    v = reshape(s(:,:,ge.k+1:ge.k+ge.B*ge.Dx),picNum,nSamp,ge.B,ge.Dx);
    for i = 1:plotnum
        secondPlotted = firstPlotted(1,i) + 1;
        if secondPlotted > ge.Dv
            secondPlotted = firstPlotted(1,i) - 1;
        end
        hor = v(:,:,:,firstPlotted(1,i));
        ver = v(:,:,:,secondPlotted);
        subplot(4,plotnum/4,i);
        scatter(hor(:),ver(:));
    end
end