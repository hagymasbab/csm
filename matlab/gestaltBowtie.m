function gestaltBowtie()
    picNum = 9;
    nSamp = 10;
    delta_x = 1;
    ge = gestaltCreate('jaguar','Dx',1024,'B',10,'obsVar',delta_x,'filters','eye','N',picNum);
    s = zeros(picNum,nSamp,ge.k+ge.B*ge.Dx);
    zs = zeros(picNum,nSamp);
    for p=1:picNum
        ge.X(p,:,:) = gestaltImageStimulus(p,B,delta_x);
    end