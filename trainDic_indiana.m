clc;
clear ;
close all;

load IndianaPines;
[r,c] = size(pixels);
tensorIndiana = reshape(pixels',ceil(sqrt(c)),ceil(sqrt(c)),r);

numExtr = 300;           % NO. of atoms in the blocks for every pic
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

sliding=8;
for i=1:r
    [blocks,idx] = my_im2col(tensorIndiana(:,:,i),[blokSize,blokSize],sliding);
    blocksExtr = randperm( length(idx),numExtr );
    blocksCluster( : , j:j+numExtr-1) = blocks( : ,blocksExtr);
    j = j+numExtr;
end
tic;
for i=1:r
    ddd=5;
    t=toc;
end

fprintf(num2str(t));
 [output] = ImageKSVD(blocksCluster,K,'blockSize',blokSize,'numKSVDIters',numIterOfKsvd,...
     'NO_atoms',L,'ModeOMP',errorFlag,'FirstAtomPre',preserveDCAtom,'reduceDC',reduceDC);
trainDic_Indiana=output.D;
save Dic_indiana_reduceMeans trainDic_Indiana_reduceMeans;