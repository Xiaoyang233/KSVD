%% initial
clear ;
close all;
clc;
load fastfirst_8_40_Dic.mat;    %data-compression rate: Dic_8L-8,Dic_4L-16,Dict_2L-32 

prompt = {'the number of atoms perserved when compression'};
lineno = 1;
defans = {'8'};
nameTitle = 'input the parameter';
answer = inputdlg(prompt,nameTitle,lineno,defans);
L = str2double(answer);

blokSize=[8,8];
K=512;

%% select Pic
%==========================================
% select pic by yourself

[PicName,PathName] = uigetfile('*.*','select the test pic');
PicStr = strcat(PathName,PicName);
%==========================================


testPicRGB = imread(PicStr);                                                            %read the test pic
if length(size(testPicRGB))==3
    testPic = rgb2gray(testPicRGB);                                                    %RGB mode to Gray mode 
    testPic = double(testPic);
else
    testPic = double(testPicRGB);
end

blocks = im2col(testPic,blokSize,'distinct');                                     %divide the pic into blocks for KSVD
blocks =double(blocks);
%% KSVD reconstruct
%-----------------------------------------------------devide line-------------------------------------------------------%
%     KSVD reconstruct

Dksvd = normc(Dksvd);
coeffice=OMP(Dksvd,blocks,L);                                                    %sparse coding
    %###############################################################################
    %encode and decode the sparse matrix

[row,cloumn,coe] = find(coeffice);

    %*************************************************
    %% encode and decode the sparse coeffice 

coe1 =round(coe);
index = find(coe1~=0);
index1 = find(coe1==0);
coe2 = coe1(index);
        %==============================================
        % encode and decode the sparse large coeffice
[totalNumLarge,tempNumLarge]=checkNum(coe2);  

[zippedCoeLarge, infoCoeLarge] = huffencode(coe2);
coe3 = huffdecode(zippedCoeLarge, infoCoeLarge);
coe3 = double(coe3);

    %*******************************************************
%% encode and decode the sparse cloumn and row
row1 = row(index);
cloumn1 = cloumn(index);

[zippedRow, infoRow] = huffencode(row1);
row2 = huffdecode(zippedRow, infoRow);

[zippedClo, infoClo] = huffencode(cloumn1);
cloumn2 = huffdecode(zippedClo, infoClo);
    %##################################################################################


[mm,nn]=size(coeffice);
coeffice1 = sparse(row2,cloumn2,coe3,mm,nn);
    %===============================================
reBlock = Dksvd*coeffice1;
reBlock1 = Dksvd*coeffice;
res=reBlock-reBlock1;
%============================================
% we must plus the mean for corresponding column because of minusing the
% mean before ksvd training 

% vecOfMeans = mean(blocks);
% reBlock = trainDic*coeffice+ones(size(blocks,1),1)*vecOfMeans;                                                            %reconstruct the blocks
%============================================

reconstKSVD = col2im(reBlock,blokSize,size(testPic),'distinct');                  %combine blocks into pic
PSNROutKsvd = 20*log10(255/sqrt(mean((reconstKSVD(:)-testPic(:)).^2)));
%-----------------------------------------------------devide line-------------------------------------------------------%

%%
%-----------------------------------------------------devide line-------------------------------------------------------%

totalCodeNo=infoCoeLarge.totalCodeNo+infoRow.totalCodeNo+infoClo.totalCodeNo;
bpp=totalCodeNo/(size(testPic,1)*size(testPic,2));

%% DCT reconstruct
%-----------------------------------------------------devide line-------------------------------------------------------%
%     DCT reconstruct

[IOut,output,BppDct] = ImageDCT(testPic,K,'L',L);
DctDic=output.D;
% reconstDCT = patch2im(IOut, size(testPic), blokSize, [1,1]);
reconstDCT = col2im(IOut,blokSize,size(testPic),'distinct');

 PSNROutDct = 20*log10(255/sqrt(mean((reconstDCT(:)-testPic(:)).^2)));
%-----------------------------------------------------devide line-------------------------------------------------------%

%% dislpay result
%-----------------------------------------------------devide line-------------------------------------------------------%
%=======================================================================================================================%
%     dislpay pervious image,the recovery image by dictionary and dictionary
%=======================================================================================================================%
%      KSVD part

figure;
subplot(2,3,1);
imshow(testPic,[]);
title('the pervious image');
subplot(2,3,2);
imshow(reconstKSVD,[]);
title(['the recovery image by K-SVD ',num2str(PSNROutKsvd),'db ',num2str(bpp),'bpp']);
subplot(2,3,3);
I = displayDictionaryElementsAsImage(Dksvd, floor(sqrt(K)), floor(size(Dksvd,2)/floor(sqrt(K))),8,8,0);
title('the K-SVD dictionary');

%=======================================================================================================================%
%  DCT part

subplot(2,3,4);
imshow(testPic,[]);
title('the pervious image');
subplot(2,3,5);
imshow(reconstDCT,[]);
title(['the recovery image by DCT ',num2str(PSNROutDct),'db ',num2str(BppDct),'bpp']);
subplot(2,3,6);
DctDic=output.D;
I = displayDictionaryElementsAsImage(DctDic, floor(sqrt(K)), floor(size(DctDic,2)/floor(sqrt(K))),8,8,0);
title('the DCT dictionary');

%-----------------------------------------------------devide line-------------------------------------------------------%
Ratio = 64/L;

testPic = uint8(testPic);
PicNameJEPG = [PicName(1:end-3),'jp2'];
jp2write(testPic,PicNameJEPG,'rate',1/Ratio,'bitdepth',8);

figure(2);
subplot(1,3,1);
imshow(testPic,[]);
title('the pervious image');

subplot(1,3,2);

JEPG2000Pic = jp2read(PicNameJEPG);
imshow(JEPG2000Pic,[]);

[PSNROutJEPG2000,~] = psnr(JEPG2000Pic,testPic);
title(['the recovery image by JEPG2000 ',num2str(PSNROutJEPG2000),'db ']);

subplot(1,3,3);
reconstKSVD = uint16(reconstKSVD);
imshow(reconstKSVD,[]);
title(['the recovery image by K-SVD ',num2str(PSNROutKsvd),'db ']);

%% comparcare
% load Dic_seterr4;
load Dic_8L;    %data-compression rate: Dic_8L-8,Dic_4L-16,Dict_2L-32 

blokSize=[8,8];
K=512;

%% select Pic
%==========================================
% select pic by yourself

[PicName,PathName] = uigetfile('*.*','select the test pic');
PicStr = strcat(PathName,PicName);
%==========================================


testPicRGB = imread(PicStr);                                                            %read the test pic
if length(size(testPicRGB))==3
    testPic = rgb2gray(testPicRGB);                                                    %RGB mode to Gray mode 
    testPic = double(testPic);
else
    testPic = double(testPicRGB);
end

blocks = im2col(testPic,blokSize,'distinct');                                     %divide the pic into blocks for KSVD
blocks =double(blocks);
%% KSVD reconstruct
%-----------------------------------------------------devide line-------------------------------------------------------%
%     KSVD reconstruct

trainDic = normc(trainDic);
coeffice=OMP(trainDic,blocks,L);                                                    %sparse coding
    %###############################################################################
    %encode and decode the sparse matrix

[row,cloumn,coe] = find(coeffice);

    %*************************************************
    %% encode and decode the sparse coeffice

coe1 =round(coe);
index = find(coe1~=0);
index1 = find(coe1==0);
coe2 = coe1(index);
        %==============================================
        % encode and decode the sparse large coeffice
[totalNumLarge,tempNumLarge]=checkNum(coe2);  

[zippedCoeLarge, infoCoeLarge] = huffencode(coe2);
coe3 = huffdecode(zippedCoeLarge, infoCoeLarge);
coe3 = double(coe3);

    %*******************************************************
%% encode and decode the sparse cloumn and row
row1 = row(index);
cloumn1 = cloumn(index);

[zippedRow, infoRow] = huffencode(row1);
row2 = huffdecode(zippedRow, infoRow);

[zippedClo, infoClo] = huffencode(cloumn1);
cloumn2 = huffdecode(zippedClo, infoClo);
    %##################################################################################

[mm,nn]=size(coeffice);
coeffice1 = sparse(row2,cloumn2,coe3,mm,nn);
% residual=coeffice-coeffice1;
% [z,w,v]=find(residual);
% disp(v);
    %===============================================
reBlock = trainDic*coeffice1;
reBlock1 = trainDic*coeffice;
res=reBlock-reBlock1;
%============================================
% we must plus the mean for corresponding column because of minusing the
% mean before ksvd training 

% vecOfMeans = mean(blocks);
% reBlock = trainDic*coeffice+ones(size(blocks,1),1)*vecOfMeans;                                                            %reconstruct the blocks
%============================================

reconstKSVD = col2im(reBlock,blokSize,size(testPic),'distinct');                  %combine blocks into pic
PSNROutKsvd = 20*log10(255/sqrt(mean((reconstKSVD(:)-testPic(:)).^2)));
%-----------------------------------------------------devide line-------------------------------------------------------%

%%
%-----------------------------------------------------devide line-------------------------------------------------------%

totalCodeNo=infoCoeLarge.totalCodeNo+infoRow.totalCodeNo+infoClo.totalCodeNo;
bpp=totalCodeNo/(size(testPic,1)*size(testPic,2));

%% DCT reconstruct
%-----------------------------------------------------devide line-------------------------------------------------------%
%     DCT reconstruct

[IOut,output,BppDct] = ImageDCT(testPic,K,'L',L);
DctDic=output.D;
% reconstDCT = patch2im(IOut, size(testPic), blokSize, [1,1]);
reconstDCT = col2im(IOut,blokSize,size(testPic),'distinct');

 PSNROutDct = 20*log10(255/sqrt(mean((reconstDCT(:)-testPic(:)).^2)));
%-----------------------------------------------------devide line-------------------------------------------------------%

%% dislpay result
%-----------------------------------------------------devide line-------------------------------------------------------%
%=======================================================================================================================%
%     dislpay pervious image,the recovery image by dictionary and dictionary
%=======================================================================================================================%
%      KSVD part

figure;
subplot(2,3,1);
imshow(testPic,[]);
title('the pervious image');
subplot(2,3,2);
imshow(reconstKSVD,[]);
title(['the recovery image by K-SVD ',num2str(PSNROutKsvd),'db ',num2str(bpp),'bpp']);
subplot(2,3,3);
I = displayDictionaryElementsAsImage(trainDic, floor(sqrt(K)), floor(size(trainDic,2)/floor(sqrt(K))),8,8,0);
title('the K-SVD dictionary');

%=======================================================================================================================%
%  DCT part

subplot(2,3,4);
imshow(testPic,[]);
title('the pervious image');
subplot(2,3,5);
imshow(reconstDCT,[]);
title(['the recovery image by DCT ',num2str(PSNROutDct),'db ',num2str(BppDct),'bpp']);
subplot(2,3,6);
DctDic=output.D;
I = displayDictionaryElementsAsImage(DctDic, floor(sqrt(K)), floor(size(DctDic,2)/floor(sqrt(K))),8,8,0);
title('the DCT dictionary');

%-----------------------------------------------------devide line-------------------------------------------------------%
% figure(2);
% imshow(testPic,[]);
% %title('the pervious image');
% saveas(gcf,'testPic.bmp');
% figure(3);
% imshow(reconstKSVD,[]);
% %title(['the recovery image by K-SVD ',num2str(PSNROutKsvd),'db ',num2str(bpp),'bpp']);
% saveas(gcf,'ReKSVD.bmp');
% figure(4);
% displayDictionaryElementsAsImage(trainDic, floor(sqrt(K)), floor(size(trainDic,2)/floor(sqrt(K))),8,8,0);
% %title('the KSVD dictionary');
% saveas(gcf,'KSVDdict.bmp');
% figure(5);
% imshow(reconstDCT,[]);
% %title(['the recovery image by DCT ',num2str(PSNROutDct),'db ',num2str(BppDct),'bpp']);
% figure(6);
% displayDictionaryElementsAsImage(DctDic, floor(sqrt(K)), floor(size(DctDic,2)/floor(sqrt(K))),8,8,0);
% %title('the DCT dictionary');
% saveas(gcf,'DCTdict.bmp');
%-----------------------------------------------------devide line-------------------------------------------------------%
Ratio = 64/L;

testPic = uint8(testPic);
% imwrite(Pic,fileName,'jp2','CompressionRatio',Ratio);
PicNameJEPG = [PicName(1:end-3),'jp2'];
jp2write(testPic,PicNameJEPG,'rate',1/Ratio,'bitdepth',8);

figure;
subplot(1,3,1);
imshow(testPic,[]);
title('the pervious image');

subplot(1,3,2);
% JEPG2000Pic = imread(fileName);

JEPG2000Pic = jp2read(PicNameJEPG);
imshow(JEPG2000Pic,[]);

[PSNROutJEPG2000,~] = psnr(JEPG2000Pic,testPic);
title(['the recovery image by JEPG2000 ',num2str(PSNROutJEPG2000),'db ']);

subplot(1,3,3);
reconstKSVD = uint16(reconstKSVD);
imshow(reconstKSVD,[]);
% title(['the recovery image by K-SVD ',num2str(PSNROutKsvd),'db ',num2str(bpp),'bpp']);
title(['the recovery image by K-SVD ',num2str(PSNROutKsvd),'db ']);