clear;
close all;
clc;
% load Dic_seterr4;
load Dic_indiana.mat;    %data-compression rate: Dic_8L-8,Dic_4L-16,Dict_2L-32 
load indianaPines.mat;

blokSize=[8,8];
K=512;
sigma = 25;

prompt = {'the band compressed','the number of atoms perserved when compression'};
lineno = 1;
defans = {'9','8'};
nameTitle = 'input the parameter';
answer = inputdlg(prompt,nameTitle,lineno,defans);

for cellnum=1:length(answer)
        parameter(cellnum) = str2double(answer{cellnum});
end
band = parameter(1);
L = parameter(2);

Pic = reshape(pixels(band,:),145,145);
blocks = im2col(Pic,blokSize,'distinct');                                     %divide the pic into blocks for KSVD
blocks =double(blocks);
%% KSVD reconstruct
%-----------------------------------------------------devide line-------------------------------------------------------%
%     KSVD reconstruct

trainDic_Indiana = normc(trainDic_Indiana);
coeffice=OMP(trainDic_Indiana,blocks,L);                                                    %sparse coding
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
        %==============================================
    %*******************************************************
%% encode and decode the sparse cloumn and row
row1 = row(index);
cloumn1 = cloumn(index);

[zippedRow, infoRow] = huffencode(row1);
row2 = huffdecode(zippedRow, infoRow);

[zippedClo, infoClo] = huffencode(cloumn1);
cloumn2 = huffdecode(zippedClo, infoClo);
    %##################################################################################
%%
[mm,nn]=size(coeffice);
coeffice1 = sparse(row2,cloumn2,coe3,mm,nn);
% residual=coeffice-coeffice1;
% [z,w,v]=find(residual);
% disp(v);
    %===============================================
reBlock = trainDic_Indiana*coeffice1;
reBlock1 = trainDic_Indiana*coeffice;
res=reBlock-reBlock1;
%============================================
% we must plus the mean for corresponding column because of minusing the
% mean before ksvd training 

% vecOfMeans = mean(blocks);
% reBlock = trainDic*coeffice+ones(size(blocks,1),1)*vecOfMeans;                                                            %reconstruct the blocks
%============================================

reconstKSVD = col2im(reBlock,blokSize,size(Pic),'distinct');                  %combine blocks into pic
PSNROutKsvd = 20*log10(255/sqrt(mean((reconstKSVD(:)-Pic(:)).^2)));
%-----------------------------------------------------devide line-------------------------------------------------------%

%%
%-----------------------------------------------------devide line-------------------------------------------------------%

% totalCodeNo=infoCoeLarge.totalCodeNo+infoCoeSmall.totalCodeNo+infoRow.totalCodeNo+infoClo.totalCodeNo;
% bpp=totalCodeNo/(size(testPic,1)*size(testPic,2));

totalCodeNo=infoCoeLarge.totalCodeNo+infoRow.totalCodeNo+infoClo.totalCodeNo;
bpp=totalCodeNo/(size(Pic,1)*size(Pic,2));

%% DCT reconstruct
%-----------------------------------------------------devide line-------------------------------------------------------%
%     DCT reconstruct

[IOut,output,BppDct] = ImageDCT(Pic,sigma,K,'L',L);
DctDic=output.D;
% reconstDCT = patch2im(IOut, size(testPic), blokSize, [1,1]);
reconstDCT = col2im(IOut,blokSize,size(Pic),'distinct');

 PSNROutDct = 20*log10(255/sqrt(mean((reconstDCT(:)-Pic(:)).^2)));
%-----------------------------------------------------devide line-------------------------------------------------------%

%% dislpay result
%-----------------------------------------------------devide line-------------------------------------------------------%
%=======================================================================================================================%
%     dislpay pervious image,the recovery image by dictionary and dictionary
%=======================================================================================================================%
%      KSVD part

figure(1);
subplot(2,3,1);
imshow(Pic,[]);
title('the pervious image');
subplot(2,3,2);
imshow(reconstKSVD,[]);
title(['the recovery image by K-SVD ',num2str(PSNROutKsvd),'db ',num2str(bpp),'bpp']);
subplot(2,3,3);
I = displayDictionaryElementsAsImage(trainDic_Indiana, floor(sqrt(K)), floor(size(trainDic_Indiana,2)/floor(sqrt(K))),8,8,0);
title('the K-SVD dictionary');

%=======================================================================================================================%
%  DCT part

subplot(2,3,4);
imshow(Pic,[]);
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
fileName = ['IndianaBand',num2str(band),'_',num2str(Ratio),'.jp2'];
Pic = uint16(Pic);
imwrite(Pic,fileName,'jp2','CompressionRatio',Ratio);
figure(2);
subplot(1,3,1);
imshow(Pic,[]);
title('the pervious image');

subplot(1,3,2);
JEPG2000Pic = imread(fileName);
PSNROutJEPG2000 = 20*log10(255/sqrt(mean((JEPG2000Pic(:)-Pic(:)).^2)));
imshow(JEPG2000Pic,[]);
title(['the recovery image by JEPG2000 ',num2str(PSNROutJEPG2000),'db ']);

subplot(1,3,3);
imshow(reconstKSVD,[]);
title(['the recovery image by K-SVD ',num2str(PSNROutKsvd),'db ',num2str(bpp),'bpp']);
%-----------------------------------------------------devide line-------------------------------------------------------%
% figure(2);
% imshow(Pic,[]);
% %title('the pervious image');
% saveas(gcf,'testPic.bmp');
% figure(3);
% imshow(reconstKSVD,[]);
% %title(['the recovery image by K-SVD ',num2str(PSNROutKsvd),'db ',num2str(bpp),'bpp']);
% saveas(gcf,'ReKSVD.bmp');
% figure(4);
% displayDictionaryElementsAsImage(trainDic_Indiana, floor(sqrt(K)), floor(size(trainDic_Indiana,2)/floor(sqrt(K))),8,8,0);
% %title('the KSVD dictionary');
% saveas(gcf,'KSVDdict.bmp');
% figure(5);
% imshow(reconstDCT,[]);
% %title(['the recovery image by DCT ',num2str(PSNROutDct),'db ',num2str(BppDct),'bpp']);
% figure(6);
% displayDictionaryElementsAsImage(DctDic, floor(sqrt(K)), floor(size(DctDic,2)/floor(sqrt(K))),8,8,0);
% %title('the DCT dictionary');
% saveas(gcf,'DCTdict.bmp');