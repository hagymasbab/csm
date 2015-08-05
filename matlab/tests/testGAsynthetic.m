setrandseed(1);
ge = gestaltCreate('temp','Dx',256,'k',3,'B',1,'N',100,'filters','OF','obsVar',0.5,'g_shape',1,'g_scale',0.1,'z_shape',2,'z_scale',2, ...
    'nullComponent',false,'generateComponents',true,'generateData',true,'componentShape','vertical-bars');
gestaltGradientAscent(ge,ge.X,10,6,'priorSamples',100,'verbose',2,'likeComp','batch','testLike',10,'learningRate',0.01,'template',false,'verbose',3,'synthetic',true,'randseed','leave','likeMethod','intuition');
