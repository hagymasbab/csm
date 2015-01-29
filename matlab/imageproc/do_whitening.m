function xw = do_whitening(x,wnode)

tlen = size(x,1);

xw=zeros(size(x));
for i=1:tlen
    xw(i,:) = (wnode*x(i,:)')';
end

% [imno patchSize] = size(x);
% 
% xw=zeros(size(x));
% for i=1:patchSize,
%     xw = xw + repmat(wnode(:,i)',imno,1).*x;
% end
