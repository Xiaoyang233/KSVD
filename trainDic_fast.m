function trainDic_fast
addpath('d:\m_sou_file\comprassive_imag_for_KSVD\ksvdbox10\');
addpath('d:\m_sou_file\comprassive_imag_for_KSVD\ksvdbox13\');



pathImag = uigetdir('d:\m_sou_file\杨老师的遥感数据\全色图像\quickBird\','the direction for train image');
if(pathImag == 0)
    return ;
end
dirImag = dir(pathImag);
NoImag = length(dirImag);

prompt = {'dictsize(the number of atoms in the dictionary)','sparsity'...
    ,'the number of interation','the number of extract patch per picture'};
lineno = 1;
defans = {'512','8','30','25000'};
title = 'input parameters';
answer = inputdlg(prompt,title,lineno,defans);

parameter = zeros(1,4);
for cellnum=1:length(answer)
        parameter(cellnum) = str2double(answer{cellnum});
end

numExtr = parameter(4);           % NO. of atoms in the blocks for every pic
j = 1;
blokSize = 8;

for i=3:NoImag
    stringDir = strcat(pathImag , '\' , dirImag(i).name);
    tempImag = imread(stringDir);
%     tempImag = rgb2gray(tempImag);
    sliding=8;
    [blocks,idx] = my_im2col(tempImag,[blokSize,blokSize],sliding);
    blocksExtr = randperm( length(idx),numExtr );
    blocksCluster( : , j:j+numExtr-1) = blocks( : ,blocksExtr);
    j = j+numExtr;
end

X = blocksCluster;

% dictionary dimensions
n = blokSize*blokSize;
m = parameter(1);

% number of examples
L = length(blocksCluster);

% sparsity of each example
k = parameter(2);

% initial dict
Pn=ceil(sqrt(m));
DCT=zeros(blokSize,Pn);
for i=0:1:Pn-1,
    V=cos([0:1:blokSize-1]'*i*pi/Pn);
    if i>0, V=V-mean(V); end;
    DCT(:,i+1)=V/norm(V);
end;
DCT=kron(DCT,DCT);

%% run k-svd training %%

params.data = X;
params.Tdata = k;
params.dictsize = m;
params.iternum = parameter(3);
params.memusage = 'high';
params.initdict = DCT(:,1:m);
params.exact = 1;
params.muthresh = 0.99;

[Dksvd,~,err] = ksvd(params,'','PARAMS','tr');
save trainDicPanliuzhi Dksvd;


%% show results %%
I = displayDictionaryElementsAsImage(Dksvd, floor(sqrt(m)), floor(size(Dksvd,2)/floor(sqrt(m))),8,8,0);

figure; plot(err); title('K-SVD error convergence');
xlabel('Iteration'); ylabel('RMSE');
printf('  Dictionary size: %d x %d', n, m);
printf('  Number of examples: %d', L);



