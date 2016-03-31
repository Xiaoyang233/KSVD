clear;
close all;
clc;
% load Dic_indiana_reduceMeans.mat;    %data-compression rate: Dic_8L-8,Dic_4L-16,Dict_2L-32 
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
imgSize = size(Pic);
dim = length(imgSize);

tic;
if(dim==3)
    j = 1;
    numblokClo = ceil(imgSize(1)/blokSize)*ceil(imgSize(2)/blokSize);
    blocks = zeros(prod(blokSizes),numblokClo*imgSize(3));
    for i = 1:imgSize(3)
        blocks(:,j:j + numblokClo - 1) = im2col(Pic(:,:,i),blokSizes,'distinct'); 
        j = j + numblokClo;
    end
else
    blocks = im2col(Pic,blokSizes,'distinct');                                     %divide the pic into blocks for KSVD
end
blocks = double(blocks);
CloumnMeans = mean(blocks);
blocks = blocks - repmat(mean(blocks),prod(blokSizes),1);
t = toc;
disp([' turn to block: ',num2str(t)]);


%% KSVD reconstruct

params.data = blocks;
params.memusage = 'high';
params.initdict = Dksvd;
Dksvd = normc(Dksvd);

L = 16;
params.Tdata = L;

tic;
% coeffice=OMP(Dksvd,blocks,L);                                                    %sparse coding
coeffice = SparseCode(params);  

t = toc;
disp([' Sparse Code: ',num2str(t)]);


tic;
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

codes.imgSize = imgSize;
codes.coefSize = size(coeffice);
codes.Min = Min;
codes.Max = Max;
codes.L = L;
codes.coepad = infoCoeLarge.pad;
codes.rowpad = infoRow.pad;
codes.symcoeLen = infoCoeLarge.symLen;
codes.symrowLen = infoRow.symLen;
codes.maxcoeLen = infoCoeLarge.maxcodelen;
codes.maxrowLen = infoRow.maxcodelen;

codes.huffTcoeIndex = infoCoeLarge.huffTIndex;
codes.huffTcoe = uint8(round(infoCoeLarge.huffT));
codes.huffTrowIndex = infoRow.huffTIndex;
codes.huffTrow = uint8(round(infoRow.huffT));

codes.bincoeLen = infoCoeLarge.totalCodeNo;
codes.coebin = logical(infoCoeLarge.string);
codes.binrowLen = infoRow.totalCodeNo;
codes.rowbin = logical(infoRow.string);

outfile = 'test.dae';
flag = daewrite(outfile, codes);

t = toc;
disp([' encode: ',num2str(t)]);