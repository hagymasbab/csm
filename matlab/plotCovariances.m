function plotCovariances(ge,dnum)
    clf;
    load iter;
    real = ge.cc;
    %xx = (1/ge.N)*(ge.X'*ge.X);
    G = ge.G(1:dnum,:);    
    V = reshape(ge.V(1:dnum,:,:),dnum*ge.B,ge.Dv);
    %size(V)
    k = size(real,2);
    n = size(pCC,2);
    s = size(VC,2);
    interleave = ceil(n/8);
    
    vertical = 2*k;
    horizontal = size(1:interleave:n,2) + 1;
    
    for i=1:interleave:n
        act = pCC{i};
        % estimated components at each step
        for j=1:k
            subplot(vertical,horizontal,(j-1)*horizontal+ceil(i/interleave));           
            viewImage(act{j},'magnif',false)
            title(sprintf('Component %d at step %d',j,i));
        end
    end
    % sample covariances at each step
    for i=1:interleave:s
        act = gVV{i};
        for j=1:k
            subplot(vertical,horizontal,(j-1+k)*horizontal+ceil(i/interleave));
            viewImage(act{j},'magnif',false)
            title(sprintf('Samp cov w. for comp. %d',j));
        end
    end
    % real components
    for j=1:k
        subplot(vertical,horizontal,j*horizontal);     
        viewImage(real{j},'magnif',false)
        title(sprintf('Real component %d',j));
    end
    % data covariance
    subplot(vertical,horizontal,(k+1)*horizontal);
    viewImage(cov(V),'magnif',false);
    title('Real hidden variable covariance');
end