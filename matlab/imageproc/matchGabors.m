function [gabors,permutation] = matchGabors(A,compare,randseed,plotStuff)
    
    close all;
    setrandseed(randseed);

    numFilt = size(A,2);
    imSize = sqrt(size(A,1));
    
    [~,maxes] = max(A);
    maxX = ceil(maxes / imSize);
    maxY = rem(maxes,imSize);
    maxY(maxY==0) = imSize;    
    
    lambdas = zeros(numFilt,1);
    orients = zeros(numFilt,1);
    phases = zeros(numFilt,1);
    sigmas = zeros(numFilt,1);
    gabors = zeros(imSize^2,numFilt);
        
    options = optimoptions('fmincon','Algorithm','interior-point','Display', 'off');
      
    permutation = zeros(numFilt,1);
    if exist(compare, 'file') == 2
        A_temp = A;
        load(compare);
        A_comp = A;
        A = A_temp;
        used = false(numFilt,1);        
    end
    
    %figure;
    for i = 1:numFilt
        printCounter(i,'stringVal','Filter','maxVal',numFilt);
        filter = A(:,i);
        
        %initparam = [4; randi([0 180]); rand()];
        if strcmp(compare,'grating')
            constraints_A = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];
            constraints_b = [imSize/2; -4; 180; 0; 1; 0];    
            initparam = [randi([4 imSize/2]); randi([0 180]); rand()];
            actfunc = @(x) difference(filter,x(1),x(2),x(3));  
            x_opt = fmincon(actfunc,initparam,constraints_A,constraints_b,[],[],[],[],[],options);       
            lambdas(i,1) = x_opt(1);
            orients(i,1) = x_opt(2);
            phases(i,1) = x_opt(3);
        elseif strcmp(compare,'gabor')
            constraints_A = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];
            constraints_b = [imSize/2; -4; 180; 0; 10; -2];    
            initparam = [randi([4 imSize/2]); randi([0 180]); randi([2 10])];
            actfunc = @(x) difference_gbr(filter,x(1),x(2),x(3),maxX(i),maxY(i)); 
            x_opt = fmincon(actfunc,initparam,constraints_A,constraints_b,[],[],[],[],[],options);       
            lambdas(i,1) = x_opt(1);
            orients(i,1) = x_opt(2);
            sigmas(i,1) = x_opt(3);
            act_gabor = gaborfilter(lambdas(i,1),orients(i,1),imSize,maxX(i),maxY(i),sigmas(i,1));
        elseif strcmp(compare,'gabor_twostep')
            constraints_A = [1 0; -1 0; 0 1; 0 -1];
            constraints_b = [imSize/2; -4; 180; 0];    
            initparam = [randi([4 imSize/2]); randi([0 180])];
            sigma0 = 10;
            actfunc = @(x) difference_gbr(filter,x(1),x(2),sigma0,maxX(i),maxY(i)); 
            x_opt = fmincon(actfunc,initparam,constraints_A,constraints_b,[],[],[],[],[],options);       
            lambdas(i,1) = x_opt(1);
            orients(i,1) = x_opt(2);
            
            sigmas(i,1) = sigma0;
%             constraints_A = [1;-1];
%             constraints_b = [10; -2];    
%             initparam = [randi([2 10])];
%             actfunc = @(x) difference_gbr(filter,lambdas(i,1),orients(i,1),x,maxX(i),maxY(i)); 
%             x_opt = fmincon(actfunc,initparam,constraints_A,constraints_b,[],[],[],[],[],options);       
%             sigmas(i,1) = x_opt;
            
            act_gabor = gaborfilter(lambdas(i,1),orients(i,1),imSize,maxX(i),maxY(i),sigmas(i,1));
            
        elseif strcmp(compare,'gabor_threestep')
            sigma0 = 10;
            lambda0 = 5;
            lambda1 = 8;
            actfunc = @(x) difference_gbr(filter,lambda0,x,sigma0,maxX(i),maxY(i)); 
            [or1,fval1] = fminbnd(actfunc,0,180);
            actfunc = @(x) difference_gbr(filter,lambda1,x,sigma0,maxX(i),maxY(i)); 
            [or2,fval2] = fminbnd(actfunc,0,180);
            if fval1 < fval2
                orients(i,1) = or1;
            else
                orients(i,1) = or2;
            end
            actfunc = @(x) difference_gbr(filter,x,orients(i,1),sigma0,maxX(i),maxY(i)); 
            lambdas(i,1) = fminbnd(actfunc,4,imSize/2);
            actfunc = @(x) difference_gbr(filter,lambdas(i,1),orients(i,1),x,maxX(i),maxY(i)); 
            sigmas(i,1) = fminbnd(actfunc,3,10);            
            %lambdas(i,1) = lambda0;
            %sigmas(i,1) = sigma0;
            act_gabor = gaborfilter(lambdas(i,1),orients(i,1),imSize,maxX(i),maxY(i),sigmas(i,1));
        else
            maxProd = -Inf;
            maxInd = 0;
            for j = 1:numFilt
                if ~used(j)
                    actProd = filter' * A_comp(:,j);
                    if actProd > maxProd
                        maxProd = actProd;
                        maxInd = j;
                    end
                end
            end
            %used(maxInd) = true;
            act_gabor = A_comp(:,maxInd);
            permutation(i) = maxInd;
        end               
        gabors(:,i) = reshape(act_gabor,imSize^2,1);
    end
    
    save('filtermatching.mat','orients','lambdas','maxX','maxY');
    
    if plotStuff
        figure;
        hist(orients,20);
        title('Orients');
        figure;
        hist(lambdas-2,20);
        title('Wavelengths');
        figure;
        hist(sigmas,20);
        title('Sigmas');
        figure;
        viewImageSet(A(:,1:100)');
        figure;
        viewImageSet(gabors(:,1:100)');

        figure;
        scatter(maxX,maxY);
        hold on;
        scatter(maxX(orients<135),maxY(orients<135));
        hold on;
        scatter(maxX(orients<90),maxY(orients<90));
        hold on;
        scatter(maxX(orients<45),maxY(orients<45));
    end
end

function rms = difference(filter,lambda,theta,phase)
    gr = grating(lambda,theta,phase,sqrt(size(filter,1)));
    rms = sum((gr(:) - filter).^2);
end

function rms = difference_gbr(filter,lambda,theta,sigma,xc,yc)
    gr = gaborfilter(lambda,theta,sqrt(size(filter,1)),xc,yc,sigma);
    rms = sum((gr(:) - filter).^2);
end

function img = gaborfilter(lambda,theta,imSize,xc,yc,sigma) 
    thetaRad = (theta / 360) * 2*pi;
    lambda = lambda - 1;
    px = xc/imSize;
    py = yc/imSize;
    [~,~,img] = gabor('theta',thetaRad,'lambda',lambda,'width',imSize,'height',imSize,'px',px,'py',py,'Sigma',sigma);             
end