function xw = do_whitening(x,wnode)

tlen = size(x,1);

xw=zeros(size(x));
for i=1:tlen
    xw(i,:) = (wnode*x(i,:)')';
end
