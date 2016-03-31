clear;
close all;
clc;
load trainDicPan16_16;
addpath('d:\m_sou_file\comprassive_imag_for_KSVD\ksvdbox10\');
blokSize = 16;
blokSizes=[blokSize,blokSize];
K=size(Dksvd,2);

[PicName,PathName,~] = uigetfile('d:\m_sou_file\杨老师的遥感数据\*.*','Selcet a file');
if(PicName == 0)
    return;
end
PicStr = [PathName,PicName];
[~, fname, ext] = fileparts(PicStr);
Pic = imread(PicStr);
if(length(size(Pic))==3)
    Pic = rgb2gray(Pic);
end

blocks = im2col(Pic,blokSizes,'distinct');                                     %divide the pic into blocks for KSVD
blocks = double(blocks);
CloumnMeans = mean(blocks);
blocks = blocks - repmat(mean(blocks),prod(blokSizes),1);


params.data = blocks;
params.memusage = 'high';
params.initdict = Dksvd;
Dksvd = normc(Dksvd);

L = 16;
params.Tdata = L;
                                                  %sparse coding
coeffice = SparseCode(params);  

[row,cloumn,coe] = find(coeffice);


[mm,nn]=size(coeffice);
coeffice1 = sparse(row,cloumn,coe,mm,nn);

reBlock = Dksvd*coeffice1;
reBlock = reBlock + repmat(CloumnMeans,prod(blokSizes),1);

reconstKSVD = col2im(reBlock,blokSizes,size(Pic),'distinct');                  %combine blocks into pic
reconstKSVD = uint8(reconstKSVD);
[PSNROutKsvd,~] = psnr(reconstKSVD,Pic);

figure(3);
imshow(reconstKSVD,[]);
title([' K-SVD (compression ratio: 16) ',num2str(PSNROutKsvd),'db ']);


L = 8;
params.Tdata = L;
                                                    %sparse coding
coeffice = SparseCode(params);  

[row,cloumn,coe] = find(coeffice);

[mm,nn]=size(coeffice);
coeffice1 = sparse(row,cloumn,coe,mm,nn);

reBlock = Dksvd*coeffice1;
reBlock = reBlock + repmat(CloumnMeans,prod(blokSizes),1);

reconstKSVD = col2im(reBlock,blokSizes,size(Pic),'distinct');                  %combine blocks into pic
reconstKSVD = uint8(reconstKSVD);
[PSNROutKsvd,~] = psnr(reconstKSVD,Pic);

figure(5);
imshow(reconstKSVD,[]);
title([' K-SVD (compression ratio: 32) ',num2str(PSNROutKsvd),'db ']);


L = 4;
params.Tdata = L;
                                                 %sparse coding
coeffice = SparseCode(params);  

[row,cloumn,coe] = find(coeffice);

[mm,nn]=size(coeffice);
coeffice1 = sparse(row,cloumn,coe,mm,nn);

reBlock = Dksvd*coeffice1;
reBlock = reBlock + repmat(CloumnMeans,prod(blokSizes),1);

reconstKSVD = col2im(reBlock,blokSizes,size(Pic),'distinct');                  %combine blocks into pic
reconstKSVD = uint8(reconstKSVD);
[PSNROutKsvd,~] = psnr(reconstKSVD,Pic);

figure(7);
imshow(reconstKSVD,[]);
title([' K-SVD (compression ratio: 64) ',num2str(PSNROutKsvd),'db ']);