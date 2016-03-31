clc;
clear all;
close all;

pathImag = uigetdir('D:\m_sou_file','the direction for train image');
dirImag = dir(pathImag);
NoImag = length(dirImag);
numExtr = 50;           % NO. of atoms in the blocks for every pic
j = 1;

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


%=======================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure;
% for i=4:NoImag
%     stringDir = strcat(pathImag , '\' , dirImag(i).name);
%     tempImag = imread(stringDir);
%     tempImag = rgb2gray(tempImag);
%     subplot(1,2,1);
%     imshow(tempImag);
%     blokSize=[16,16];
%     sliding=16;
%     [blocks,idx] = my_im2col(tempImag,blokSize,sliding);
%     
%     num = 0;
%     temp1Imag = zeros(size(tempImag));
%      [rows,cols] = ind2sub((size(tempImag)-blokSize+1),idx);     %get the index of first pixal in the block
%          for numCol = 1:sliding:length(cols)                                 %look the variable table of clos,you will see it
%             for numRow = 1:blokSize(1)                                         %look the variable table of rows,you will see it
%            num = num+1;                                                             %the value of the num will set 0 everytime you read the next pic 
%            tempBlok=reshape(blocks(:,num), blokSize);
%             temp1Imag( rows(numRow) : rows(numRow)+blokSize(1)-1 , cols(numCol) : cols(numCol)+blokSize(2)-1 ) = tempBlok;
%             end
%          end
%     subplot(1,2,2);
%     imshow(temp1Imag,[]);                                                          %the type of  entry's value in "temp1Imag" is double,so we must add "[]" into imshow
%     pause(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i=4:NoImag
    stringDir = strcat(pathImag , '\' , dirImag(i).name);
    tempImag = imread(stringDir);
    tempImag = rgb2gray(tempImag);
    sliding=8;
    [blocks,idx] = my_im2col(tempImag,[blokSize,blokSize],sliding);
    blocksExtr = randperm( length(idx),numExtr );
    blocksCluster( : , j:j+numExtr-1) = blocks( : ,blocksExtr);
    j = j+numExtr;
end
 [output] = ImageKSVD(blocksCluster,K,'blockSize',blokSize,'numKSVDIters',numIterOfKsvd,...
     'NO_atoms',L,'ModeOMP',errorFlag,'FirstAtomPre',preserveDCAtom,'reduceDC',reduceDC);
trainDic=output.D;
save Dic trainDic;