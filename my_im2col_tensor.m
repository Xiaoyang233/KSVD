function [blocks,idx] = my_im2col_tensor(I,blkSize,slidingDis)
% description: arrange 3 dimension image to block
% input:
% I: image
% blkSize: 2 dimension vector which imply the size of bolck
% slidingDis: the slideing size in the image
% 
% output:
% blocks: the arrange image 
% idx: the first index in each blocks

[r,c,d] = size(I);
idxMat = zeros([r,c]-blkSize+1);
idxMat([[1:slidingDis:end-1],end],[[1:slidingDis:end-1],end]) = 1; % take blocks in distances of 'slidingDis', but always take the first and last one (in each row and column).
idx = find(idxMat);
[rows,cols] = ind2sub(size(idxMat),idx);
blocks = zeros(prod(blkSize),length(idx),d);
for i = 1:length(idx)
    currBlock = I(rows(i):rows(i)+blkSize(1)-1,cols(i):cols(i)+blkSize(2)-1,:);
    blocks(:,i,:) = reshape(currBlock,prod(blkSize),d);
end
