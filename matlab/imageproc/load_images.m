function rawImages=load_images

w=1536;h=1024;

d0=pwd;

cd /Users/karaj/csnl/vanhateren

% load -ascii filenames.dat

fid=fopen('filenames.dat','r');
A=fread(fid,inf, 'uint8=>char');
fclose(fid);

A=reshape(A,13, round(length(A')/13))';
A=A(:,1:12);

rawImages=zeros(size(A,1),w,h);
% rawImages=zeros(w,h,size(A,1));
% fprintf('Loading images...\n')
% fprintf('%3d',1);
for iid=1:size(A,1)
%     fprintf('\b\b\b%3d',iid);
    f1=fopen(A(iid,:),'rb','ieee-be');
    rawImages(iid,:,:)=permute(fread(f1,[w,h],'uint16'),[3 1 2]);
%     rawImages(:,:,iid)=fread(f1,[w,h],'uint16');
    fclose(f1);
end
% fprintf('\n');

cd(d0);

% rawImages=permute(rawImages,[3 2 1]);
