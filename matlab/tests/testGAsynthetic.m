function testGAsynthetic(Dx,k,learningRate,N_all,batchSize,priorSamples,testLike,noiseLevel)

    setrandseed(1);    
    ge = gestaltCreate('temp','Dx',Dx,'k',k,'B',1,'N',N_all,'filters','OF','obsVar',noiseLevel,'g_shape',1,'g_scale',0.1,'z_shape',2,'z_scale',2, ...
        'nullComponent',false,'generateComponents',true,'generateData',true,'componentShape','vertical-bars');
    
    gestaltGradientAscent(ge,ge.X,batchSize,8,'priorSamples',priorSamples,'verbose',2,'likeComp','batch','testLike',testLike, ...
        'learningRate',learningRate,'template',false,'verbose',2,'synthetic',true,'randseed','leave','likeMethod','intuition','initSigma',noiseLevel);
end