setrandseed(1);
ge = gestaltCreate('temp','Dx',256,'k',3,'B',1,'N',10000,'filters','OF','obsVar',0.5,'g_shape',1,'g_scale',0.1,'z_shape',2,'z_scale',2, ...
    'nullComponent',false,'generateComponents',true,'generateData',true,'componentShape','vertical-bars');
gestaltGradientAscent(ge,ge.X,800,8,'priorSamples',200,'verbose',2,'likeComp','batch','testLike',200,'learningRate',0.01,'template',false,'verbose',2,'synthetic',true,'randseed','leave','likeMethod','intuition');
