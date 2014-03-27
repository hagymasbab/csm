function plotCovariances(ge,dnum)
    clf;
    load iter;
    real = ge.cc;
    %xx = (1/ge.N)*(ge.X'*ge.X);
    G = ge.G(1:dnum,:);
    V = ge.V(1:dnum,:);
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
            imshow(act{j})
            title(sprintf('Component %d at step %d',j,i));
        end
    end
    % sample covariances at each step
    for i=1:interleave:s
        subplot(k+1,len,k*len+ceil(i/interleave));
        imshow(VC{i})
        title(sprintf('Sample covariance at step %d',i));
    end
    % real components
    for j=1:k
        subplot(k+1,len,j*len);     
        imshow(real{j})
        title(sprintf('Real component %d',j));
    end
    % data covariance
    subplot(k+1,len,(k+1)*len);
    imshow(cov(V));
    title('Real hidden variable covariance');
    
    
end