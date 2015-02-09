function gabors = matchGabors(A,compare,randseed)
    
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
    gabors = zeros(numFilt);
        
    options = optimoptions('fmincon','Algorithm','interior-point','Display', 'off');
        
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
            
            constraints_A = [1;-1];
            constraints_b = [10; -2];    
            initparam = [randi([2 10])];
            actfunc = @(x) difference_gbr(filter,lambdas(i,1),orients(i,1),x,maxX(i),maxY(i)); 
            x_opt = fmincon(actfunc,initparam,constraints_A,constraints_b,[],[],[],[],[],options);       
            sigmas(i,1) = x_opt;
            
            act_gabor = gaborfilter(lambdas(i,1),orients(i,1),imSize,maxX(i),maxY(i),sigmas(i,1));
        end               
        gabors(:,i) = reshape(act_gabor,imSize^2,1);
    end
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

function rms = difference(filter,lambda,theta,phase)
    gr = grating(lambda,theta,phase,sqrt(size(filter,1)));
    rms = sum((gr(:) - filter).^2);
end

function rms = difference_gbr(filter,lambda,theta,sigma,xc,yc)
    gr = gaborfilter(lambda,theta,sqrt(size(filter,1)),xc,yc,sigma);
    rms = sum((gr(:) - filter).^2);
end

function img = grating(lambda,theta,phase,imSize)
    % imSize/2 >= lambda >= 3
    % 180 > theta >= 0
    % 1 > phase >= 0
    freq = imSize/lambda;
    phaseRad = (phase * 2* pi);
    X0 = 1:imSize;
    X0 = (X0 / imSize) - .5;
    [Xm,Ym] = meshgrid(X0, X0); 
    thetaRad = (theta / 360) * 2*pi;        % convert theta (orientation) to radians
    Xt = Xm * cos(thetaRad);                % compute proportion of Xm for given orientation
    Yt = Ym * sin(thetaRad);                % compute proportion of Ym for given orientation
    XYt = Xt + Yt;                          % sum X and Y components
    XYf = XYt * freq * 2*pi;                % convert to radians and scale by frequency
    img = sin( XYf + phaseRad);    
end

function img = gaborfilter(lambda,theta,imSize,xc,yc,sigma) 
    thetaRad = (theta / 360) * 2*pi;
    lambda = lambda - 1;
    px = xc/imSize;
    py = yc/imSize;
    [~,~,img] = gabor('theta',thetaRad,'lambda',lambda,'width',imSize,'height',imSize,'px',px,'py',py,'Sigma',sigma);             
end