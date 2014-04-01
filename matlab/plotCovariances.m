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
    len = size(1:interleave:n,2) + 1;
    for i=1:interleave:n
        act = pCC{i};
        % estimated components at each step
        for j=1:k
            subplot(k+1,len,(j-1)*len+ceil(i/interleave));           
            viewImage(act{j},'magnif',false)
            title(sprintf('Component %d at step %d',j,i));
        end
    end
    % sample covariances at each step
    for i=1:interleave:s
        subplot(k+1,len,k*len+ceil(i/interleave));
        viewImage(VC{i},'magnif',false)
        title(sprintf('Sample covariance at step %d',i));
    end
    % real components
    for j=1:k
        subplot(k+1,len,j*len);     
        viewImage(real{j},'magnif',false)
        title(sprintf('Real component %d',j));
    end
    % data covariance
    subplot(k+1,len,(k+1)*len);
    viewImage(cov(V),'magnif',false);
    title('Real hidden variable covariance');
end