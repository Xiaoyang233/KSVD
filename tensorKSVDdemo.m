
clc
clear;

pathImag = uigetdir('D:\m_sou_file','the direction for train image');
if pathImag==0
	return
end
dirImag = dir(pathImag);
NoImag = length(dirImag);
numExtr = 50;           % NO. of atoms in the blocks for every pic
j = 0;

%======================================================
%input of prarmeter
% prompt = {'the number of KSVD interation','the size of blocks(input only one number)',...
%     'the number of atoms in KSVD interation','the mode of OMP(1 is set error mode)',...
%     'perserve the first atom(1 to perserve)','minus means of blocks(1 to minus)'};
% lineno = 1;
% defans = {'25','8','8','0','1','0'};
% title = 'input the parameter of KSVD';
% answer = inputdlg(prompt,title,lineno,defans);
% for cellnum=1:length(answer)
%         parameter(cellnum) = str2double(answer{cellnum});
% end
% %=======================================================

% K=512;                      %dictionary size 
% L=parameter(3);
% numIterOfKsvd = parameter(1);
% blokSize=parameter(2);
% preserveDCAtom = parameter(5);
% errorFlag =  parameter(4);
% reduceDC = parameter(6);

blockSize = 8;
epsilon = 0.02;
% numExtr = 50;
numIteration = 20;
selectNum = 2;

for i=4:NoImag
    stringDir = strcat(pathImag , '\' , dirImag(i).name);
    tempImag = imread(stringDir);
    % sliding=blockSize;
    % [blocks,idx] = my_im2col_tensor(tempImag,[blockSize,blockSize],sliding);
    % blocksExtr = randperm( length(idx),numExtr );
    % blocksCluster( : , j:j+numExtr-1,:) = blocks( : ,blocksExtr,:);
    % j = j+numExtr;
    tempBlock = randomBlock(tempImag,[blockSize,blockSize],selectNum);
    for k =1:selectNum
    	blocksCluster(:,j*selectNum*blockSize + 1 + (k-1)*blockSize:j*selectNum*blockSize + 1 + k*blockSize - 1,:) = tempBlock{k};
    end
    j = j + 1;
end
colNum = 50;
rowNum = 40;
blocksCluster = reshape(blocksCluster,rowNum*blockSize,colNum*blockSize,3);
blocksCluster = double(blocksCluster);

Pad = zeros(4,colNum*blockSize,3);
blocksCluster = [blocksCluster;Pad];

blockImgSize = size(blocksCluster);
D = cell(1,3);
Sparisty = [0.016, 0.016, 0.09];
DictSize = blockImgSize*8;
S = round(DictSize.*Sparisty);
% Create an initial dictionary from the DCT frame
for i = 1:3
	Pn = ceil(sqrt(DictSize(i)));
	Xn = ceil(sqrt(blockImgSize(i)));
	D{i}=zeros(Xn,Pn);
	for k=0:1:Pn-1,
	    V=cos([0:1:Xn-1]'*k*pi/Pn);
	    if k>0, V=V-mean(V); end;
	    D{i}(:,k+1)=V/norm(V);
	end;
	D{i}=kron(D{i},D{i});
end

D{3} = D{3}(1:3,:);

[X,D] = tensorKSVD1(D,blocksCluster,S,epsilon,numIteration);
