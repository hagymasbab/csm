function secondOrder(imgpath)
    close all;
    
    im1 = zeros(64,64);
    im1(31:54,31:34) = 2;
    im1(30,32:33) = 2;
    im1(55,32:33) = 2;
    viewImage(im1)
    print('-dpng',[imgpath '/comp1.png']);
    
    im2 = im1';
    figure;
    viewImage(im2)
    print('-dpng',[imgpath '/comp2.png']);
    
    im3 = fliplr(im2);
    figure;
    viewImage(im3)
    print('-dpng',[imgpath '/comp3.png']);
    
    figure;
    imc = im1+im2+im3;
    viewImage(imc,'useMax',true);    
    hold on
    scatter(32,32,'rs','filled')
    print('-dpng',[imgpath '/comp_c.png']);
end