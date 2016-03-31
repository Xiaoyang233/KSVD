clear;
close all;
clc;
% load Dic_indiana_reduceMeans.mat;    %data-compression rate: Dic_8L-8,Dic_4L-16,Dict_2L-32 
load trainDicPan;

addpath('d:\m_sou_file\comprassive_imag_for_KSVD\ksvdbox10\');
blokSize = 8;
blokSizes=[blokSize,blokSize];
K=size(Dksvd,2);
% addpath('d:\m_sou_file\comprassive_imag_for_KSVD\ksvdbox10\');

% prompt = {'the band compressed','the number of atoms perserved when compression'};
% lineno = 1;
% defans = {'9','8'};
% nameTitle = 'input the parameter';
% answer = inputdlg(prompt,nameTitle,lineno,defans);
% 
% parameter = zeros(2,1);
% for cellnum=1:length(answer)
%         parameter(cellnum) = str2double(answer{cellnum});
% end
% band = parameter(1);
% L = parameter(2);

[PicName,PathName,~] = uigetfile('d:\m_sou_file\杨老师的遥感数据\*.*','Selcet a file');
if(PicName == 0)
    return;
end
PicStr = [PathName,PicName];
[~, fname, ext] = fileparts(PicStr);
Pic = imread(PicStr);

f = figure();
scrsz = get(0,'ScreenSize');
imshow(Pic,'border','tight','initialmagnification','fit');
set(f,'Position',[0 30 scrsz(3) scrsz(4)-95]);

[y,x] = ginput();
close;

Pic = subim(Pic, 256, 256, ceil(x), ceil(y));
Pic = uint8(Pic);

if(length(size(Pic))==3)
    Pic = rgb2gary(Pic);
end

blocks = im2col(Pic,blokSizes,'distinct');                                     %divide the pic into blocks for KSVD
blocks = double(blocks);
CloumnMeans = mean(blocks);
blocks = blocks - repmat(mean(blocks),prod(blokSizes),1);



%% KSVD reconstruct

params.data = blocks;
params.memusage = 'high';
params.initdict = Dksvd;
Dksvd = normc(Dksvd);

L = 8;
params.Tdata = L;


% coeffice=OMP(Dksvd,blocks,L);                                                    %sparse coding
coeffice = SparseCode(params);  

tic;
[row,cloumn,coe] = find(coeffice);

coe1 =round(coe);
index = find(coe1~=0);
coe2 = coe1(index);

[zippedCoeLarge, infoCoeLarge] = huffencode(coe2);
coe3 = huffdecode(zippedCoeLarge, infoCoeLarge);
coe3 = double(coe3);

row1 = row(index);
cloumn1 = cloumn(index);

[zippedRow, infoRow] = huffencode(row1);
row2 = huffdecode(zippedRow, infoRow);

[zippedClo, infoClo] = huffencode(cloumn1);

cloumn2 = huffdecode(zippedClo, infoClo);
t(1) = toc;

tic;
% coeffice =round(coeffice);
% [r_encode,infor] = matrix2huff(coeffice);
% coeffice2 = huff2matrix(r_encode,infor);
% coeffice2 = sparse(coeffice2);

[~,~,coe] = find(coeffice);

coe1 =round(coe);
index = find(coe1~=0);
coe21 = coe1(index);

Min = min(coe21(:));
Max = max(coe21(:));
coe2 = linerGary(coe21,0,255)+1;

[zippedCoeLarge, infoCoeLarge] = huffencode(coe2);
coe4 = huffdecode(zippedCoeLarge, infoCoeLarge);

coe4 = double(coe4);
coe4 = linerGary(coe4-1,Min,Max);
coe4 = round(coe4);

disp(isequal(coe21(:),coe4(:)));

t(2) = toc;

[mm,nn]=size(coeffice);
coeffice1 = sparse(row2,cloumn2,coe3,mm,nn);

reBlock = Dksvd*coeffice1;
reBlock = reBlock + repmat(CloumnMeans,prod(blokSizes),1);

coeffice2 = sparse(row2,cloumn2,coe4,mm,nn);
reBlock1 = Dksvd*coeffice2;
reBlock1 = reBlock1 + repmat(CloumnMeans,prod(blokSizes),1);
%============================================
% we must plus the mean for corresponding column because of minusing the
% mean before ksvd training 

% vecOfMeans = mean(blocks);
% reBlock = trainDic*coeffice+ones(size(blocks,1),1)*vecOfMeans;                                                            %reconstruct the blocks
%============================================

reconstKSVD = col2im(reBlock,blokSizes,size(Pic),'distinct');                  %combine blocks into pic
reconstKSVD = uint8(reconstKSVD);
[PSNROutKsvd,~] = psnr(reconstKSVD,Pic);

reconstKSVD1 = col2im(reBlock1,blokSizes,size(Pic),'distinct');                  %combine blocks into pic
reconstKSVD1 = uint8(reconstKSVD1);
[PSNROutKsvd1,~] = psnr(reconstKSVD1,Pic);

%-----------------------------------------------------devide line-------------------------------------------------------%

% totalCodeNo=infoCoeLarge.totalCodeNo+infoCoeSmall.totalCodeNo+infoRow.totalCodeNo+infoClo.totalCodeNo;
% bpp=totalCodeNo/(size(testPic,1)*size(testPic,2));

totalCodeNo=infoCoeLarge.totalCodeNo+infoRow.totalCodeNo+infoClo.totalCodeNo;
bpp=totalCodeNo/(size(Pic,1)*size(Pic,2));



%-----------------------------------------------------devide line-------------------------------------------------------%
Ratio = 64/L;
fileName = [fname,'_',num2str(Ratio),'.jp2'];
jp2write(Pic,fileName,'rate',1/Ratio,'bitdepth',8);

figure(1);
subplot(1,4,1);
imshow(Pic,[]);
title('the pervious image');

subplot(1,4,2);
JEPG2000Pic = jp2read(fileName);
imshow(JEPG2000Pic,[]);

[PSNROutJEPG2000,~] = psnr(JEPG2000Pic,Pic);
title([' JEPG2000 (compression ratio: 8) ',num2str(PSNROutJEPG2000),'db ']);

subplot(1,4,3);

imshow(reconstKSVD,[]);
% title(['the recovery image by K-SVD ',num2str(PSNROutKsvd),'db ',num2str(bpp),'bpp']);
title([' K-SVD (compression ratio: 8) ',num2str(PSNROutKsvd),'db ']);

subplot(1,4,4);

imshow(reconstKSVD1,[]);
title([' K-SVD (compression ratio: 8) ',num2str(PSNROutKsvd1),'db ']);

% =================================================================================================================

L = 4;

params.Tdata = L;

% coeffice=OMP(Dksvd,blocks,L);                                                    %sparse coding
coeffice = SparseCode(params);                                                      %sparse coding

[row,cloumn,coe] = find(coeffice);

coe1 =round(coe);
index = find(coe1~=0);
coe2 = coe1(index);

tic;
[zippedCoeLarge, infoCoeLarge] = huffencode(coe2);
coe3 = huffdecode(zippedCoeLarge, infoCoeLarge);
coe3 = double(coe3);

row1 = row(index);
cloumn1 = cloumn(index);

[zippedRow, infoRow] = huffencode(row1);
row2 = huffdecode(zippedRow, infoRow);

[zippedClo, infoClo] = huffencode(cloumn1);
cloumn2 = huffdecode(zippedClo, infoClo);

[mm,nn]=size(coeffice);
coeffice1 = sparse(row2,cloumn2,coe3,mm,nn);
t(3) = toc;

tic;
% [r_encode,infor] = matrix2huff(coeffice);
% coeffice2 = huff2matrix(r_encode,infor);
% coeffice2 = sparse(coeffice2);

[~,~,coe] = find(coeffice);

coe1 =round(coe);
index = find(coe1~=0);
coe2 = coe1(index);

Min = min(coe2(:));
Max = max(coe2(:));
coe2 = linerGary(coe2,1,256);

[zippedCoeLarge, infoCoeLarge] = huffencode(coe2);
coe4 = huffdecode(zippedCoeLarge, infoCoeLarge);
coe4 = double(coe4);
coe4 = linerGary(coe4,Min,Max);

t(4) = toc;

reBlock = Dksvd*coeffice1;
reBlock = reBlock + repmat(CloumnMeans,prod(blokSizes),1);

coeffice2 = sparse(row2,cloumn2,coe4,mm,nn);
reBlock1 = Dksvd*coeffice2;
reBlock1 = reBlock1 + repmat(CloumnMeans,prod(blokSizes),1);

reconstKSVD = col2im(reBlock,blokSizes,size(Pic),'distinct');                  %combine blocks into pic
reconstKSVD = uint8(reconstKSVD);
[PSNROutKsvd,~] = psnr(reconstKSVD,Pic);

reconstKSVD1 = col2im(reBlock1,blokSizes,size(Pic),'distinct');                  %combine blocks into pic
reconstKSVD1 = uint8(reconstKSVD1);
[PSNROutKsvd1,~] = psnr(reconstKSVD1,Pic);

%-----------------------------------------------------devide line-------------------------------------------------------%

% totalCodeNo=infoCoeLarge.totalCodeNo+infoCoeSmall.totalCodeNo+infoRow.totalCodeNo+infoClo.totalCodeNo;
% bpp=totalCodeNo/(size(testPic,1)*size(testPic,2));

totalCodeNo=infoCoeLarge.totalCodeNo+infoRow.totalCodeNo+infoClo.totalCodeNo;
bpp=totalCodeNo/(size(Pic,1)*size(Pic,2));



%-----------------------------------------------------devide line-------------------------------------------------------%
Ratio = 64/L;
fileName = [fname,'_',num2str(Ratio),'.jp2'];
jp2write(Pic,fileName,'rate',1/Ratio,'bitdepth',8);

figure(2);
subplot(1,4,1);
imshow(Pic,[]);
title('the pervious image');

subplot(1,4,2);
JEPG2000Pic = jp2read(fileName);
imshow(JEPG2000Pic,[]);

[PSNROutJEPG2000,~] = psnr(JEPG2000Pic,Pic);
title([' JEPG2000 (compression ratio: 16)',num2str(PSNROutJEPG2000),'db ']);

subplot(1,4,3);

imshow(reconstKSVD,[]);
% title(['the recovery image by K-SVD ',num2str(PSNROutKsvd),'db ',num2str(bpp),'bpp']);
title([' K-SVD (compression ratio: 16)',num2str(PSNROutKsvd),'db ']);

subplot(1,4,4);

imshow(reconstKSVD1,[]);
title([' K-SVD (compression ratio: 8) ',num2str(PSNROutKsvd1),'db ']);

% =================================================================================================================

L = 2;

params.Tdata = L;

% coeffice=OMP(Dksvd,blocks,L);                                                    %sparse coding
coeffice = SparseCode(params);                                                      %sparse coding

tic;
[row,cloumn,coe] = find(coeffice);

coe1 =round(coe);
index = find(coe1~=0);
index1 = find(coe1==0);
coe2 = coe1(index);

[zippedCoeLarge, infoCoeLarge] = huffencode(coe2);
coe3 = huffdecode(zippedCoeLarge, infoCoeLarge);
coe3 = double(coe3);

row1 = row(index);
cloumn1 = cloumn(index);

[zippedRow, infoRow] = huffencode(row1);
row2 = huffdecode(zippedRow, infoRow);

[zippedClo, infoClo] = huffencode(cloumn1);
cloumn2 = huffdecode(zippedClo, infoClo);

[mm,nn]=size(coeffice);
coeffice1 = sparse(row2,cloumn2,coe3,mm,nn);

t(5) = toc;

tic;
% [r_encode,infor] = matrix2huff(coeffice);
% coeffice2 = huff2matrix(r_encode,infor);
% coeffice2 = sparse(coeffice2);

[~,~,coe] = find(coeffice);

coe1 =round(coe);
index = find(coe1~=0);
coe2 = coe1(index);

Min = min(coe2(:));
Max = max(coe2(:));
coe2 = linerGary(coe2,1,256);

[zippedCoeLarge, infoCoeLarge] = huffencode(coe2);
coe4 = huffdecode(zippedCoeLarge, infoCoeLarge);
coe4 = double(coe4);
coe4 = linerGary(coe4,Min,Max);

t(6) = toc;

reBlock = Dksvd*coeffice1;
reBlock = reBlock + repmat(CloumnMeans,prod(blokSizes),1);

coeffice2 = sparse(row2,cloumn2,coe4,mm,nn);
reBlock1 = Dksvd*coeffice2;
reBlock1 = reBlock1 + repmat(CloumnMeans,prod(blokSizes),1);

reconstKSVD = col2im(reBlock,blokSizes,size(Pic),'distinct');                  %combine blocks into pic
reconstKSVD = uint8(reconstKSVD);
[PSNROutKsvd,~] = psnr(reconstKSVD,Pic);

reconstKSVD1 = col2im(reBlock1,blokSizes,size(Pic),'distinct');                  %combine blocks into pic
reconstKSVD1 = uint8(reconstKSVD1);
[PSNROutKsvd1,~] = psnr(reconstKSVD1,Pic);
%-----------------------------------------------------devide line-------------------------------------------------------%

% totalCodeNo=infoCoeLarge.totalCodeNo+infoCoeSmall.totalCodeNo+infoRow.totalCodeNo+infoClo.totalCodeNo;
% bpp=totalCodeNo/(size(testPic,1)*size(testPic,2));

totalCodeNo=infoCoeLarge.totalCodeNo+infoRow.totalCodeNo+infoClo.totalCodeNo;
bpp=totalCodeNo/(size(Pic,1)*size(Pic,2));



%-----------------------------------------------------devide line-------------------------------------------------------%
Ratio = 64/L;
fileName = [fname,'_',num2str(Ratio),'.jp2'];
jp2write(Pic,fileName,'rate',1/Ratio,'bitdepth',8);

figure(3);
subplot(1,4,1);
imshow(Pic,[]);
title('the pervious image');

subplot(1,4,2);
JEPG2000Pic = jp2read(fileName);
imshow(JEPG2000Pic,[]);

[PSNROutJEPG2000,~] = psnr(JEPG2000Pic,Pic);
title([' JEPG2000 (compression ratio: 32) ',num2str(PSNROutJEPG2000),'db ']);

subplot(1,4,3);

imshow(reconstKSVD,[]);
% title(['the recovery image by K-SVD ',num2str(PSNROutKsvd),'db ',num2str(bpp),'bpp']);
title([' K-SVD (compression ratio: 32) ',num2str(PSNROutKsvd),'db ']);

subplot(1,4,4);

imshow(reconstKSVD1,[]);
title([' K-SVD (compression ratio: 8) ',num2str(PSNROutKsvd1),'db ']);

figure(4);
I = displayDictionaryElementsAsImage(Dksvd, floor(sqrt(K)), floor(size(Dksvd,2)/floor(sqrt(K))),blokSize,blokSize,0);