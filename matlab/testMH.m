function [samples,rr] = testMH(dim,nSamp,sampler,varargin)
    p = inputParser;
    addParamValue(p,'sphere',false,@islogical);
    addParamValue(p,'verbose',0,@isnumeric);
    addParamValue(p,'thin',0,@isnumeric);
    addParamValue(p,'burnin',0,@isnumeric);
    addParamValue(p,'propVar',0.1,@isnumeric);
    addParamValue(p,'adapt',false,@islogical);
    addParamValue(p,'covariance',zeros(dim,dim));
    addParamValue(p,'metnum',0,@isnumeric);
    addParamValue(p,'wide',true,@islogical);
    addParamValue(p,'transform',false,@islogical);
    parse(p,varargin{:});
    burnin = p.Results.burnin;
    thin = p.Results.thin;
    
    fprintf('Creating covariance matrix\n');
    if max(max(p.Results.covariance)) > 0
        C = p.Results.covariance;
    elseif p.Results.sphere
        C = 2 * eye(dim);
    elseif p.Results.wide       
        wide = false;
        while ~wide
            b = rand(dim,dim);
            C = 0.1*eye(dim) + (b'*b)/dim;
            wide = true;            
            [~,S,~] = svd(C);
            if max(diag(S)) / min(diag(S)) > 1e+3
                wide = false;
            end           
        end
    else
        b = rand(dim,dim);
        C = (b'*b)/dim;
    end
    iC = inv(C);
    mu = -3*ones(dim,1);
    %logpdf   = @(x) (-1/2) * (dim*log(2*pi) + log(det(C)) + (x-mu)' * iC * (x-mu) );
    logpdf   = @(x) (-1/2) * ( log(det(C)) + (x-mu)' * iC * (x-mu) );
    grad = @(x) -iC * (x - mu);
    init = 2 * rand(dim,1);
    propCov = p.Results.propVar * eye(dim);
    fprintf('Starting sampler\n');
    if strcmp(sampler,'mh')
        [samples,rr] = metropolisHastings(init,logpdf,propCov,nSamp,burnin,thin,'verbose',p.Results.verbose,'adapt',p.Results.adapt);
    elseif strcmp(sampler,'hmc')
        subplot(2,1,1)
        plotGaussian(mu,C);
        hold on
        
        [~,S,~] = svd(C);
        leapsteps = ceil(max(diag(S)) / min(diag(S)))
        stepsize = 0.9*min(diag(S))
        %leapsteps = 10000;
        %stepsize = 0.006;
        
        [samples,rr] = hamiltonianMC(init,logpdf,grad,nSamp,leapsteps,stepsize,'verbose',p.Results.verbose,'plot',true);
        %figure
    elseif strcmp(sampler,'cmc')
        m = p.Results.metnum;
        dims = [];
        if m>0
            dims = 1:m;
        end
        
        plot = false;
        if m>1 || dim-m>1
            plot = true;
            clf;
            subplot(1,2,1);
            if m > 1
                plotGaussian(mu(1:m,:),C(1:m,1:m))
                title('Metropolis');
            end
            hold on
            subplot(1,2,2);
            if dim-m>1
                plotGaussian(mu(m+1:dim,:),C(m+1:dim,m+1:dim))
                title('Hamilton');
            end
            hold on
        end 
        modgrad = @(x) grad(x).* [zeros(m,1);ones(dim-m,1)];
        bounds = [];
        %bounds = [1 -3 2;2 -3 2];
        [samples,rr] = combinedMC(init,dims,logpdf,modgrad,nSamp,C(m+1:dim,m+1:dim),p.Results.propVar,'verbose',p.Results.verbose,'plot',plot,'bounds',bounds,'transformMomentum',p.Results.transform);
    else
        fprintf('Invalid sampler specified\n');
        return;
    end
    
%     figure
%     subplot(2,1,1);
%     plotGaussian(mu,C);
%     hold on;
%     scatter(samples(:,1),samples(:,2));
%     xl = xlim();
%     yl = ylim();
%     hold off;
%     subplot(2,1,2);
%     smoothhist2D(samples,5,[100, 100]);
%     xlim(xl);
%     ylim(yl);
%     set(gca,'YDir','normal')
    
    function plotGaussian(mu,C)
        [x,y] = meshgrid(-5:0.1:3,-5:0.1:3);
        pr = mvnpdf([x(:) y(:)], mu(1:2,:)', C(1:2,1:2));
        contour(x,y,reshape(pr,size(x,1),size(x,1)));
    end
end