clear;
close all;
clc;

% ×Öµä 256*1024
load trainDicPan16_16;
blokSize = 16;
blokSizes=[blokSize,blokSize];
K=size(Dksvd,2);

[PicName,PathName,~] = uigetfile('.\*.*','Selcet a file');
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

Dksvd = normc(Dksvd);
% 16±¶Ñ¹Ëõ
L = 16;
coeffice=OMP(Dksvd,blocks,L);                                                    %sparse coding 

[row,cloumn,coe] = find(coeffice);

coe1 =round(coe);
index = find(coe1~=0);
coe2 = coe1(index);

% ÏßÐÔ±ä»» 0-255£¨Á¿»¯£© 
Min = min(coe2(:));
Max = max(coe2(:));
coe21 = linerGary(coe2,0,255)+1;

%huffman±àÂë
[zippedCoeLarge, infoCoeLarge] = huffencode(coe21);

row1 = row(index);
[zippedRow, infoRow] = huffencode(row1);

cloumn1 = cloumn(index);

% huffman½âÂë
coe3 = huffdecode(zippedCoeLarge, infoCoeLarge);
coe3 = double(coe3);
coe3 = linerGary(coe3-1,Min,Max);
coe3 = round(coe3);

row1 = huffdecode(zippedRow, infoRow);

disp(isequal(coe2(:),coe3(:)));


[mm,nn]=size(coeffice);
coeffice1 = sparse(row1,cloumn1,coe3,mm,nn);

reBlock = Dksvd*coeffice1;
reBlock = reBlock + repmat(CloumnMeans,prod(blokSizes),1);

reconstKSVD = col2im(reBlock,blokSizes,size(Pic),'distinct');                  %combine blocks into pic
reconstKSVD = uint8(reconstKSVD);
[PSNROutKsvd,~] = psnr(reconstKSVD,Pic);

figure(1);
imshow(Pic,[]);
title('the pervious image');

figure(2);
imshow(reconstKSVD,[]);
title([' K-SVD (compression ratio: 16) ',num2str(PSNROutKsvd),'db ']);


% 32±¶Ñ¹Ëõ
L = 8;
coeffice=OMP(Dksvd,blocks,L);                                                    %sparse coding 

[row,cloumn,coe] = find(coeffice);

coe1 =round(coe);
index = find(coe1~=0);
coe2 = coe1(index);

Min = min(coe2(:));
Max = max(coe2(:));
coe21 = linerGary(coe2,0,255)+1;

[zippedCoeLarge, infoCoeLarge] = huffencode(coe21);

row1 = row(index);
[zippedRow, infoRow] = huffencode(row1);

cloumn1 = cloumn(index);

coe3 = huffdecode(zippedCoeLarge, infoCoeLarge);
coe3 = double(coe3);
coe3 = linerGary(coe3-1,Min,Max);
coe3 = round(coe3);

row1 = huffdecode(zippedRow, infoRow);

disp(isequal(coe2(:),coe3(:)));


[mm,nn]=size(coeffice);
coeffice1 = sparse(row1,cloumn1,coe3,mm,nn);

reBlock = Dksvd*coeffice1;
reBlock = reBlock + repmat(CloumnMeans,prod(blokSizes),1);

reconstKSVD = col2im(reBlock,blokSizes,size(Pic),'distinct');                  %combine blocks into pic
reconstKSVD = uint8(reconstKSVD);
[PSNROutKsvd,~] = psnr(reconstKSVD,Pic);

figure(3);
imshow(reconstKSVD,[]);
title([' K-SVD (compression ratio: 32) ',num2str(PSNROutKsvd),'db ']);

% 64±¶Ñ¹Ëõ
L = 4;

coeffice=OMP(Dksvd,blocks,L);                                                    %sparse coding 

[row,cloumn,coe] = find(coeffice);

coe1 =round(coe);
index = find(coe1~=0);
coe2 = coe1(index);

Min = min(coe2(:));
Max = max(coe2(:));
coe21 = linerGary(coe2,0,255)+1;

[zippedCoeLarge, infoCoeLarge] = huffencode(coe21);

row1 = row(index);
[zippedRow, infoRow] = huffencode(row1);

cloumn1 = cloumn(index);


coe3 = huffdecode(zippedCoeLarge, infoCoeLarge);
coe3 = double(coe3);
coe3 = linerGary(coe3-1,Min,Max);
coe3 = round(coe3);

row1 = huffdecode(zippedRow, infoRow);

disp(isequal(coe2(:),coe3(:)));

[mm,nn]=size(coeffice);
coeffice1 = sparse(row1,cloumn1,coe3,mm,nn);

reBlock = Dksvd*coeffice1;
reBlock = reBlock + repmat(CloumnMeans,prod(blokSizes),1);

reconstKSVD = col2im(reBlock,blokSizes,size(Pic),'distinct');                  %combine blocks into pic
reconstKSVD = uint8(reconstKSVD);
[PSNROutKsvd,~] = psnr(reconstKSVD,Pic);

figure(4);
imshow(reconstKSVD,[]);
title([' K-SVD (compression ratio: 64) ',num2str(PSNROutKsvd),'db ']);

