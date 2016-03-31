clc;
clear all;
close all;

% pathImag = uigetdir('D:\m_sou_file','the direction for train image');
% dirImag = dir(pathImag);
% NoImag = length(dirImag);
numExtr = 361;           % NO. of atoms in the blocks for every pic
j = 1;
load IndianaPines;

%======================================================
%input of prarmeter
prompt = {'the number of KSVD interation','the size of blocks(input only one number)',...
    'the number of atoms in KSVD interation','the mode of OMP(1 is set error mode)',...
    'perserve the first atom(1 to perserve)','minus means of blocks(1 to minus)'};
lineno = 1;
defans = {'25','8','8','0','1','0'};
title = 'input the parameter of KSVD';
answer = inputdlg(prompt,title,lineno,defans);
for cellnum=1:length(answer)
        parameter(cellnum) = str2double(answer{cellnum});
end
%=======================================================

K=512;                      %dictionary size 
L=parameter(3);
numIterOfKsvd = parameter(1);
blokSize=parameter(2);
preserveDCAtom = parameter(5);
errorFlag =  parameter(4);
reduceDC = parameter(6);
sigma=25;

%=======================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for i=4:NoImag
% % parfor i=4:NoImag
%     stringDir = strcat(pathImag , '\' , dirImag(i).name);
%     tempImag = imread(stringDir);
%     tempImag = rgb2gray(tempImag);
% %     subplot(1,2,1);
% %     imshow(tempImag);
% 
%     sliding=8;
%     [blocks,idx] = my_im2col(tempImag,[blokSize,blokSize],sliding);
%     blocksExtr = randperm( length(idx),numExtr );
%     blocksCluster( : , j:j+numExtr-1) = blocks( : ,blocksExtr);
%     j = j+numExtr;
% end
dim=min(size(pixels));
for i=1:dim
    temp=reshape(pixels(i,:),[145,145]);
    blocks = im2col(temp,[blokSize,blokSize],'distinct');
    blocksCluster( : , j:j+numExtr-1) = blocks;
    j = j+numExtr;
end
[output] = ImageKSVD(blocksCluster,sigma,K,'blockSize',blokSize,'numKSVDIters',numIterOfKsvd,...
     'NO_atoms',L,'ModeOMP',errorFlag,'FirstAtomPre',preserveDCAtom);
trainDic=output.D;
save Dic trainDicMulSpec;